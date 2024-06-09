"""A script for comparing the statistical power of the Qini and AUTOC.

This script contains code for comparing the Qini and AUTOC on a level
playing field. Because the Qini Coefficient weights do not necessarily
sum to 1 while those of the AUTOC always do, it can be difficult to
compare the power of these two metrics without centering and scaling
the Qini Coefficient weights. We do that here for a class of natural,
simple functions governing the potential outcomes given treatment
and covariates.
"""

import argparse
import os
from typing import Callable, Iterable, Mapping, Optional, Tuple, List

import econml
from econml import grf
import matplotlib as mpl
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sklearn
from sklearn.model_selection import KFold
from tqdm import tqdm


def scaled_RATE(sorted_scores: np.ndarray, method: str = "AUTOC") -> float:
    """Calculate a centered, scaled estimate of the RATE

    Calculates a centered and scaled version of the Qini (and AUTOC).
    Given some fixed sample size, n, the variance of the weights
    applied to each doubly robust score in the AUTOC is constrained
    to be 1 while the variance of the weights applied to these same
    scores for the Qini coefficient tends toward 0.5 for large n. In this
    function, we rescale the Qini coefficient weights to also have variance
    1 in order to fairly compare the two metrics. While this rescaling
    changes the value of the point estimate for these methods, it does
    not change the statistical power of the two approaches.

    Args:
        sorted_scores: A [num_samples, ] array containing "scores" (nearly
            unbiased proxies of the CATE, such as IPW or AIPW scores), sorted
            in order of the priority scores output by the prioritization rule
            under consideration.
        method: A string indicating which of the two methods ("AUTOC" vs.
            "QINI") should be calculated

    Returns:
        RATE (float): A point estimate of the Rank-Weighted Average Treatment
            Effect, the exact form of which is specified by `method`.
    """
    n = len(sorted_scores)
    H: np.ndarray = np.array([])
    if method == "AUTOC":
        rev_inv_rank = 1.0 / np.arange(n, 0, -1)
        # AUTOC should be centered/scaled by default
        H = np.flip(np.cumsum(rev_inv_rank)) - 1
    elif method == "QINI":
        rev_sort_u = np.arange(n, 0, -1) / n
        # We explicitly center and scale the Qini coefficient
        H = (rev_sort_u - np.mean(rev_sort_u)) / np.std(rev_sort_u)
    else:
        raise ValueError(f"Method {method} is not supported")
    RATE = np.mean(H * sorted_scores)
    return RATE


def aipw_func(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    e: Optional[float] = 0.5,
    m: Optional[Callable] = None,
    n_folds: int = 2,
) -> np.ndarray:
    """Estimates Augmented Inverse Propensity Weight (AIPW) scores

    Args:
        X: A [num_samples, num_features] numpy array with shape representing
            N samples and p input covariates
        Y: A [num_samples, 1] numpy array representing outcomes corresponding
            to N individuals given by X.
        W: A [num_samples, 1] numpy array with W[i] = 1 if the ith individual
            was treated and 0 otherwise. Note that only binary treatments are
            currently supported.
        e: A float that represents the probability the ith subject would be
            treated (i.e., prob that W[i] = 1). We assume known propensity
            scores in this case.
        m: A function that takes in an [num_samples, num_features + 1] matrix
            representing the covariates for each subject and a final column
            representing treatment assignment and returns the estimated
            marginal outcome. If this is defined, it is used instead of
            training new nuisance parameter models/marginal outcome estimators
            from scratch. Defining m to be the ground truth marginal response
            curves can be used to calculate Oracle scores.

    Returns:
        AIPW_Scores: A [N, 1] numpy array representing the estimated
            AIPW scores for each individual
    """
    n = X.shape[0]
    e_hat = np.repeat(e, n)

    # Shape handling to make sure everything works with eg sklearn, econml
    if len(W.shape) == 1:
        W = W.reshape(-1, 1)
    if len(X.shape) == 1:
        X = X.reshape(-1, 1)
    Y = Y.flatten()

    # AIPW scores for continuous and binary outcomes
    if m is not None:  # If m is known then we generate oracle scores
        mu_hat_1 = m(X, np.ones(n))
        mu_hat_0 = m(X, np.zeros(n))
    else:
        mu_hat_1 = np.ones((n, 1)) * -np.infty
        mu_hat_0 = np.ones((n, 1)) * -np.infty
        kf = sklearn.model_selection.KFold(n_splits=n_folds, shuffle=True)
        for train_idx, test_idx in kf.split(X):
            n_test = len(test_idx)

            # Use an honest forest to estimate baseline model/outcomes
            outcome_model = econml.grf.RegressionForest()
            outcome_model.fit(np.hstack([X[train_idx], W[train_idx]]), Y[train_idx])

            # Predict outcomes under treatment for the held out individuals
            mu_hat_1[test_idx] = outcome_model.predict(
                np.hstack([X[test_idx], np.ones((n_test, 1))])
            )

            # Predict outcomes under control for the held out individuals
            mu_hat_0[test_idx] = outcome_model.predict(
                np.hstack([X[test_idx], np.zeros((n_test, 1))])
            )

        # Make sure we gave an estimate for every subject
        assert not np.any(mu_hat_1 == -np.infty)
        assert not np.any(mu_hat_0 == -np.infty)

    # Use standard formula to estimate AIPW scores for each subject
    AIPW_scores = (
        mu_hat_1.flatten()
        - mu_hat_0.flatten()
        + W.flatten() / e_hat.flatten() * (Y - mu_hat_1.flatten())
        - (1 - W.flatten()) / (1 - e_hat.flatten()) * (Y - mu_hat_0.flatten())
    )

    return AIPW_scores


def ipw_func(
    X: np.ndarray, Y: np.ndarray, W: np.ndarray, e: Optional[float] = 0.5
) -> np.ndarray:
    """Estimates Inverse Propensity Weight (IPW) scores

    Args:
        X: A [num_samples, num_features] numpy array with shape representing
            N samples each with d-dimensional input covariates
        Y: A [num_samples, 1] numpy array representing outcomes corresponding
            to N individuals with covariates given by X.
        W: A [num_samples, 1] numpy array with W[i] = 1 if the ith individual
            was treated and 0 otherwise. Note that only binary treatments are
            currently supported.
        e: A float that represents the probability the ith subject would be
            treated (i.e., prob that W[i] = 1). We assume known propensity
            scores in this case.

    Returns:
        IPW_Scores: A [N, 1] numpy array representing the estimated AIPW scores
            for each individual
    """
    e_hat = np.repeat(e, X.shape[0])
    W = W.flatten()
    Y = Y.flatten()
    e_hat = e_hat.flatten()
    IPW_scores = (W * Y / e_hat) - ((1.0 - W) * Y / (1.0 - e_hat))
    return IPW_scores


def get_scores(X, Y, W, e=0.5, m=None, scoring_type="AIPW") -> np.ndarray:
    """Estimate scores (proxies for the CATE) for each subject

    Args:
        X: A [num_samples, num_features] numpy array with shape representing
            N samples each with d-dimensional input covariates
        Y: A [num_samples, 1] numpy array representing outcomes corresponding
            to N individuals with covariates given by X.
        W: A [num_samples, 1] numpy array with W[i] = 1 if the ith individual
            was treated and 0 otherwise. Note that only binary treatments are
            currently supported.
        e: A float that represents the probability the ith subject would be
            treated (i.e., prob that W[i] = 1)
        m: A function that takes in an [num_samples, num_features + 1] matrix
            representing the covariates for each subject and a final column
            representing treatment assignment and returns the estimated
            marginal outcome. If this is defined, it is used instead of
            training new nuisance parameter models/marginal outcome estimators
            from scratch. Defining m to be the ground truth marginal response
            curves can be used to calculate Oracle scores.
        scoring_type: One of 'IPW', 'AIPW', or 'Oracle', representing inverse
            propensity weighted scores, augmented inverse propensity weighted
            scores, or 'Oracle' scores which are AIPW scores where the true
            conditional response surfaces are known and given. Note that
            the 'Oracle' scores are not the exact individualized treamtent
            difference, but rather the true expected difference in potential
            outcomes conditioned on covariates.

    Returns:
        A [num_samples, ] numpy array where the ith element of the array
        represents either an IPW, AIPW, or Oracle score for that subject.
        These scores are nearly unbiased proxies for the CATE.
    """
    if scoring_type == "IPW":
        return ipw_func(X, Y, W, e)
    elif scoring_type == "AIPW":
        return aipw_func(X, Y, W, e)
    elif scoring_type == "Oracle" and m is not None:
        return aipw_func(X, Y, W, e, m)
    elif scoring_type == "Oracle" and m is None:
        raise ValueError(
            "Cannot provide oracle scores without ground "
            "truth marginal response function!"
        )
    else:
        raise ValueError(f"Scoring method {scoring_type} not supported")


def response_fn(X: np.ndarray, W: np.ndarray, t: float) -> np.ndarray:
    """Generate treatment- and covariate-dependent response

    This function generates outcomes conditioned on covariates and
    treatment for each sample in 1, ..., num_samples. Specifically,
    if treatment W[i] = 0 then outcome m[i] = 0. If W[i] = 1, then
    m[i] = max(-2X[i] / t^2 + 2/t, 0). This gives a response surface
    under treatment for which subjects with X < t have a nonzero
    difference in expected outcome for control vs. treatment and
    subjects with X ≥ t have zero difference (i.e., CATE = 0).
    The scaling makes it so that the integral under the curve
    for treated individuals is the same irrespective of which
    value of t is chosen.

    Args:
        X: A [num_samples, 1] numpy array with shape representing
            N samples, each with 1-dimensional input covariates
        W: A [num_samples, 1] numpy array with W[i] = 1 if the ith individual
            was treated and 0 otherwise. Note that only binary treatments are
            currently supported.
        t: A threshold parameter for the response surface. Potential outcomes
            under control always have expectation 0 independent of covariates,
            in this particular example; however, subjects with covariate X < t
            will have potential outcome under treatment Y(1) be non-zero, while
            subjects with covariate X ≥ t will have potential outcome Y(1) be
            zero.

    Returns:
        m: A numpy array where the ith element of the array represents the
            expected difference in potential outcomes for the ith subject
            conditioned on that subject's covariates and actual treatment
            assignemnt W[i].
    """
    m = np.empty(X.shape[0])
    if not np.all(W == 1):
        m[W == 0] = 0
    if not np.all(W == 0):
        m[W == 1] = np.maximum(-2.0 / t**2.0 * X.flatten() + 2.0 / t, 0)[W == 1]
    return m


def generate_data(
    N: int, p: float = 0.5, frac_zero_cate: float = 0.5, seed: int = 42
) -> Tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, Callable
]:
    """Generate synthetic covariates and potential outcomes for each subject

    This function synthesizes samples with 1-dimensional covariates for X
    and treatment assignment W = 1 (treated) or W = 0 (control) assigned
    randomly to each subject with probability p. Expected potential outcomes
    under treatment and control are generated via the `response_fn`.
    Potential outcomes Y0 (for potential outcome under control) and
    Y1 (for potential outcome under treatment) are generated via the
    `response_fn` but with i.i.d. gaussian noise for each subject.
    The expected difference in potential outcomes conditioned on covariates
    (or "CATE") is calculated as the difference in the potential outcomes
    before adding i.i.d. noise.

    Args:
        N : int
            The number of samples to generate
        p : float
            The probability of treatment for each sample.
        frac_zero_cate : float
            The proportion of samples for whom the difference in potential
            outcomes is exactly zero (see `response_fn`)
        seed : int
            Random seed to enable reproducibility

    Returns:
        X: A numpy array of shape [N, ] representing covariates for
            each subject
        W: Numpy array of shape [N, ] where W[i] = 1 if ith simulated subject
            is treated, 0 otherwise
        Y: Numpy array of shape [N, ] representing the observed outcomes
            for each simulated subject
        Y0: Numpy array of shape [N, ] representing the potential outcome
            under control (no treatment) for each subject
        Y1: Numpy array of shape [N, ] representing potential outcome
            under treatment for each subject
        tau: Numpy array of shape [N, ] representing expected difference
            in potential outcomes conditioned on covariates, for each
            subject's covariates
         m_fn: A callable (python function) that takes in arrays `x` and `w`
            and returns the expected outcome conditioned on covariates/
            treatment assignment.
    """
    np.random.seed(seed)
    X = np.random.uniform(size=N)
    W = (np.random.uniform(size=N) < p) * 1
    m_fn = lambda X, W: response_fn(X, W, t=frac_zero_cate)
    m0 = m_fn(X, np.zeros(X.shape[0]))
    m1 = m_fn(X, np.ones(X.shape[0]))
    Y0 = m0 + np.random.randn(N) * 0.2
    Y1 = m1 + np.random.randn(N) * 0.2
    Y = np.empty(N)
    Y[:] = np.NaN
    Y[W == 0] = Y0[W == 0]
    Y[W == 1] = Y1[W == 1]
    tau = m1 - m0
    return X, W, Y, Y0, Y1, tau, m_fn


def plot_score_comparisons(
    qini_estimate_vecs: List[np.ndarray],
    autoc_estimate_vecs: List[np.ndarray],
    score_names: List[str],
    save_dir: str,
    fname: str,
) -> None:
    """
    Generates and saves a plot comparing the distribution of the
    estimated Rank-weighted Average Treatment Effect (RATE)
    for different scoring types (e.g., IPW, AIPW, and Oracle scores).

    Args:
        qini_estimate_vecs (List[np.ndarray]): List of arrays containing QINI
            coefficient estimates for different scoring types.
        autoc_estimate_vecs (List[np.ndarray]): List of arrays containing AUTOC
            estimates for different scoring types.
        score_names (List[str]): List of names for the different scoring types.
        save_dir (str): Directory to save the resulting plot.
        fname (str): Filename for the saved plot.

    Returns:
        Nothing, but saves an image in "save_dir/fname"
    """
    score_color_map = {"IPW": "tab:blue", "AIPW": "tab:orange", "Oracle": "tab:green"}
    # Initialize plot settings to produce LaTeX-style text
    mpl.rcParams["axes.formatter.use_mathtext"] = True
    mpl.rcParams["font.family"] = "serif"
    cmfont = font_manager.FontProperties(
        fname=mpl.get_data_path() + "/fonts/ttf/cmr10.ttf"
    )
    mpl.rcParams["font.serif"] = cmfont.get_name()
    mpl.rcParams["mathtext.fontset"] = "cm"
    mpl.rcParams["axes.unicode_minus"] = False

    # Set up the figure and axis
    fig, axes = plt.subplots(1, 2, figsize=(5, 3), dpi=300)

    for qini_estimates, autoc_estimates, score_name in zip(
        qini_estimate_vecs, autoc_estimate_vecs, score_names
    ):
        # Left plot
        sns.histplot(
            autoc_estimates,
            ax=axes[0],
            label=score_name,
            kde=True,
            stat="density",
            color=score_color_map[score_name],
            edgecolor="none",
            bins=50,
        )

        # Right plot
        sns.histplot(
            qini_estimates,
            ax=axes[1],
            label=score_name,
            kde=True,
            stat="density",
            color=score_color_map[score_name],
            edgecolor="none",
            bins=50,
        )

    # Left plot
    axes[0].set_title("AUTOC")
    axes[0].set_xlabel("RATE")
    axes[0].set_ylabel("Density")

    # Right plot
    axes[1].set_title("QINI Coefficient")
    axes[1].set_xlabel("RATE")
    axes[1].set_ylabel("")

    # Legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=1,
        bbox_to_anchor=(0.5, 0.95),
        title="Scoring",
        facecolor="white",
        fontsize=10,
        framealpha=1,
        edgecolor="black",
    )

    plt.tight_layout()
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(os.path.join(save_dir, fname))
    plt.close("all")


def plot_qini_vs_autoc(
    qini_autoc_pairs: List[Tuple[np.ndarray, np.ndarray]],
    X: np.ndarray,
    tau: np.ndarray,
    Y0: np.ndarray,
    Y1: np.ndarray,
    save_dir: str,
    fname: str,
    scoring_types: List[str] = None,
    combine_plots: bool = False,
) -> None:
    """Plot AUTOC and Qini comparison together with potential outcomes given X

    Args:
        qini_autoc_pairs : List[Tuple[np.ndarray, np.ndarray]]
            Each element in the list corresponds to a separate scoring function
            (e.g., IPW scores, Oracle scores, AIPW scores).
            For each pair of arrays in the list, the first array is
            Qini estimates from B experiments, the second is AUTOC estimates
            from B experiments, where B is the number of simulations performed.

        X : np.ndarray
            N samples, each with 1-D input covariates.

        tau : np.ndarray
            Expected difference in potential outcomes for each sample

        Y0 : np.ndarray
            Potential outcomes under control (no treatment) for each subject

        Y1 : np.ndarray
            Potential outcomes under treatment for each subject

        save_dir : str
            The directory in which the plot figures should be saved

        fname : str
            The filename that should be given to the saved plot

        scoring_types : List[str]
            The scoring types

        combine_plots : bool
            If True, combine all plots into a single figure with four rows

    Returns:
        Nothing, but saves an image in "save_dir/fname"
    """
    # Initialize plot settings to produce LaTeX-style text
    mpl.rcParams["axes.formatter.use_mathtext"] = True
    mpl.rcParams["font.family"] = "serif"
    cmfont = font_manager.FontProperties(
        fname=mpl.get_data_path() + "/fonts/ttf/cmr10.ttf"
    )
    mpl.rcParams["font.serif"] = cmfont.get_name()
    mpl.rcParams["mathtext.fontset"] = "cm"
    mpl.rcParams["axes.unicode_minus"] = False

    # Define the number of rows for the plots, 2 rows for histograms + 1 row for scatter plot
    num_rows = len(qini_autoc_pairs) + 1 if combine_plots else 2

    fig, axs = plt.subplots(num_rows, figsize=(3, 5 * num_rows // 2), dpi=300)

    min_RATE = np.infty
    max_RATE = -np.infty
    for qini_estimates, autoc_estimates in qini_autoc_pairs:
        tmp_min = min(np.min(qini_estimates), np.min(autoc_estimates))
        if tmp_min < min_RATE:
            min_RATE = tmp_min
        tmp_max = max(np.max(qini_estimates), np.max(autoc_estimates))
        if tmp_max > max_RATE:
            max_RATE = tmp_max
    hist_bins = np.linspace(min_RATE, max_RATE, num=50)

    for idx, ((qini_estimates, autoc_estimates), scoring_type) in enumerate(
        zip(qini_autoc_pairs, scoring_types)
    ):
        ax_hist = axs[idx] if combine_plots else axs[0]

        sns.histplot(
            qini_estimates,
            kde=True,
            stat="density",
            label="Qini Coefficient (Linear)",
            ax=ax_hist,
            color="tab:blue",
            edgecolor="none",
            bins=hist_bins,
        )
        sns.histplot(
            autoc_estimates,
            kde=True,
            stat="density",
            label="AUTOC (Logarithmic)",
            ax=ax_hist,
            color="tab:orange",
            edgecolor="none",
            bins=hist_bins,
        )

        ax_hist.legend(
            title="Weighting",
            bbox_to_anchor=(0.5, 1.4),
            loc="upper center",
            edgecolor="black",
            facecolor="white",
            framealpha=1.0,
            fontsize=10,
        )
        x_label = r"$\widehat{RATE}$ Distribution"
        if scoring_type:
            x_label += f" ({scoring_type} scores)"
        ax_hist.set_xlabel(x_label, fontsize=10)
        ax_hist.set_ylabel("Density", fontsize=10)
        ax_hist.set_xlim(min_RATE, max_RATE)

    # Plot scatter plot only once, either as the last subplot or separately if not combined
    ax_scatter = axs[-1]
    sorted_indices = np.argsort(X)
    ax_scatter.scatter(
        X[sorted_indices],
        (Y1 - Y0)[sorted_indices],
        alpha=1.0,
        s=15,
        linewidth=0.5,
        facecolors="none",
        edgecolors="#999999",
        label=r"$Y(1) - Y(0)$",
    )
    ax_scatter.plot(
        X[sorted_indices],
        tau[sorted_indices],
        linewidth=1.0,
        color="#984ea3",
        alpha=0.8,
        label=r"$\tau(X)$",
    )
    # ax_scatter.set_title("Individual Treatment Effects and CATE")
    ax_scatter.set_xlabel("$X$", fontsize=10)
    ax_scatter.set_ylabel("Treatment Effect", fontsize=10)
    ax_scatter.legend(edgecolor="black", facecolor="white", fontsize=10)

    # fig.subplots_adjust(hspace=0.1)
    plt.tight_layout()
    os.makedirs(save_dir, exist_ok=True)
    plt.savefig(os.path.join(save_dir, fname))
    plt.close("all")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n_sims",
        type=int,
        default=10000,  # TODO: Reset to 10k
        help="Number of simulations to run "
        "(each simulation results in one RATE estimate).",
    )
    parser.add_argument(
        "--n_samples",
        type=int,
        default=400,
        help="Number of samples to generate within each simulation.",
    )
    parser.add_argument(
        "--p_treat",
        type=float,
        default=0.5,
        help="Proportion of samples that are assigned to treatment arm.",
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default="figures/",
        help="Directory in which generated images should be saved",
    )
    parser.add_argument(
        "--frac_zero_cate",
        nargs="+",
        default=[1.0, 0.5, 0.1],
        help='List where each element "l" represents an experiment '
        'in which "l" out of 1000 subjects have CATE = 0, on avg',
    )
    parser.add_argument(
        "--weighting_functions",
        nargs="+",
        choices=["QINI", "AUTOC"],
        default=["QINI", "AUTOC"],
        help="List of weighting functions to use. Acceptable values are 'QINI' and 'AUTOC'.",
    )
    parser.add_argument(
        "--scoring_types",
        nargs="+",
        choices=["IPW", "AIPW", "Oracle"],
        default=["IPW", "AIPW", "Oracle"],
        help="List of scoring types to use. Acceptable values are 'IPW', 'AIPW', and 'Oracle'.",
    )
    args = parser.parse_args()

    # Initialize the results dictionary
    RATE_estimates = {}
    for frac in args.frac_zero_cate:
        for weight_fn in args.weighting_functions:
            for scoring_type in args.scoring_types:
                RATE_estimates[(frac, weight_fn, scoring_type)] = []

    # Simulations take about 75 minutes on a MacBook Pro w/
    # a 2.6 GHz 6-core Intel Core i7 processor and 16 GB 2667 MHz DDR4 RAM
    for frac in args.frac_zero_cate:
        print(f"Running simulations with {frac * 100:.0f}% non-zero CATE")
        for b in tqdm(range(args.n_sims)):
            X, W, Y, Y0, Y1, tau, m_fn = generate_data(
                N=args.n_samples,
                p=args.p_treat,
                frac_zero_cate=frac,
                seed=b,
            )
            for scoring_type in args.scoring_types:
                scores = get_scores(
                    X, Y, W, e=args.p_treat, m=m_fn, scoring_type=scoring_type
                )
                for weight_fn in args.weighting_functions:
                    RATE_estimates[(frac, weight_fn, scoring_type)] += [
                        scaled_RATE(scores[np.argsort(X)], method=weight_fn)
                    ]

        # Generate a reference set of data with a particular seed
        X, W, _, Y0, Y1, tau, m_fn = generate_data(
            N=args.n_samples, p=args.p_treat, frac_zero_cate=frac, seed=0
        )

        # Generate the pairs for each scoring type
        qini_autoc_pairs = [
            (
                RATE_estimates[(frac, "QINI", scoring_type)],
                RATE_estimates[(frac, "AUTOC", scoring_type)],
            )
            for scoring_type in args.scoring_types
        ]

        # Filename for the combined plot
        fname = f"qini_vs_autoc_when_{int(frac * 100)}_pct_have_nonzero_cate.png"

        # Call the modified plot function with the combined flag set to True
        plot_qini_vs_autoc(
            qini_autoc_pairs=qini_autoc_pairs,
            X=X,
            tau=tau,
            Y0=Y0,
            Y1=Y1,
            save_dir=args.out_dir,
            fname=fname,
            scoring_types=args.scoring_types,
            combine_plots=True,  # Set to True to combine all plots into a single figure
        )

        plot_score_comparisons(
            qini_estimate_vecs=[
                RATE_estimates[(frac, "QINI", score)] for score in args.scoring_types
            ],
            autoc_estimate_vecs=[
                RATE_estimates[(frac, "AUTOC", score)] for score in args.scoring_types
            ],
            score_names=args.scoring_types,
            save_dir=args.out_dir,
            fname=f"score_comparisons_for_{int(frac * 100)}_pct_nonzero.png",
        )
