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
import econml
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib as mpl
import numpy as np
import os
import seaborn as sns
from tqdm import tqdm
from typing import Callable, Iterable, Mapping, Optional, Tuple


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
            treated (i.e., prob that W[i] = 1)
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
    e_hat = np.repeat(e, X.shape[0])

    # SHAPE HANDLING
    if len(W.shape) == 1:
        W = W.reshape(-1, 1)
    if len(X.shape) == 1:
        X = X.reshape(-1, 1)
    Y = Y.flatten()

    # AIPW scores for continuous and binary outcomes
    if m is not None:
        treated_preds = m(X, np.ones(X.shape[0]))
        control_preds = m(X, np.zeros(X.shape[0]))
        actual_preds = m(X, W.flatten())
    else:
        outcome_hf = econml.grf.RegressionForest()
        outcome_hf.fit(
            np.hstack([X, W]), Y
        )  # Use an honest forest to estimate baseline model/outcomes
        treated_preds = outcome_hf.predict(
            np.hstack([X, np.ones_like(W)])
        )  # Estimate outcome under treatment
        control_preds = outcome_hf.predict(
            np.hstack([X, np.zeros_like(W)])
        )  # Estimate outcome under control
        actual_preds = outcome_hf.predict(
            np.hstack([X, W])
        )  # Estimate outcomes under true/obs treatment assignment

    # Standard formula to estimate AIPW scores
    AIPW_scores = (
        treated_preds.flatten()
        - control_preds.flatten()
        + ((W.flatten() - e_hat.flatten()) * (Y - actual_preds.flatten()))
        / (e_hat * (1 - e_hat))
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
            treated (i.e., prob that W[i] = 1)

    Returns:
        IPW_Scores: A [N, 1] numpy array representing the estimated AIPW scores
            for each individual
    """
    e_hat = np.repeat(e, X.shape[0])
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
    N: int, p: float = 0.5, seed: int = 42
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
        N: An integer representing the number of samples to generate
        p: A float repreesnting the probability of treatment for each sample.
            Also used as the proportion of samples for whom the difference in
            potential outcomes is exactly zero (see `response_fn`)
        seed: An (integer) random seed to enable reproducibility

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
    m_fn = lambda X, W: response_fn(X, W, t=p)
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


def plot_auc_results(
    auc_vec: Mapping[str, Iterable],
    X: np.ndarray,
    tau: np.ndarray,
    Y0: np.ndarray,
    Y1: np.ndarray,
    save_dir: str,
    fname: str,
) -> None:
    """Plot AUTOC and Qini comparison together with potential outcomes given X

    Args:
        auc_vec: A [num_samples, ] numpy array with shape representing
            the RATE with element order corresponding to X
        X: A [num_samples, ] numpy array with shape representing
            N samples, each with 1-dimensional input covariates
        tau: A [num_samples, ] numpy array representing the expected
            difference in potential outcomes given covariates for
            each sample
        Y0: A [num_samples, ] numpy array representing the potential
            outcome under control (no treatment) for each subject
        Y1: A [num_samples, ] numpy array representing the potential
            outcome under treatment for each subject
        save_dir: A string representing the directory in which the
            plot figures should be saved
        fname: The filename that should be given to the saved plot

    Returns:
        Nothing, but saves an image in "save_dir/fname"
    """
    # Initialize plot settings to produce LaTeX-style text
    mpl.rcParams["font.family"] = "serif"
    cmfont = font_manager.FontProperties(
        fname=mpl.get_data_path() + "/fonts/ttf/cmr10.ttf"
    )
    mpl.rcParams["font.serif"] = cmfont.get_name()
    mpl.rcParams["mathtext.fontset"] = "cm"
    mpl.rcParams["axes.unicode_minus"] = False

    weight_type_dict = {
        "QINI": "Qini Coefficient (Linear)",
        "AUTOC": "AUTOC (Logarithmic)",
    }
    fig, axs = plt.subplots(2, figsize=(3, 5), dpi=300)
    for weight_type in ["QINI", "AUTOC"]:
        sns.distplot(
            auc_vec[weight_type], label=weight_type_dict[weight_type], ax=axs[0]
        )

    axs[0].legend(
        title="Weighting",
        bbox_to_anchor=(0.5, 1.4),
        loc="upper center",
        edgecolor="black",
        facecolor="white",
        framealpha=1.0,
        fontsize=10,
    )
    axs[0].set_xlabel(r"$\widehat{RATE}$  Distribution", fontsize=10)
    axs[0].set_ylabel("Density", fontsize=10)
    axs[1].scatter(
        X[np.argsort(X)],
        (Y1 - Y0)[np.argsort(X)],
        alpha=1.0,
        s=15,
        linewidth=0.5,
        facecolors="none",
        edgecolors="#999999",
        label=r"$Y(1) - Y(0)$",
    )
    axs[1].plot(
        X[np.argsort(X)],
        tau[np.argsort(X)],
        linewidth=1.0,
        color="#984ea3",
        alpha=0.8,
        label=r"$\tau(X)$",
    )
    axs[1].set_xlabel("$X$", fontsize=10)
    axs[1].set_ylabel("Treatment Effect", fontsize=10)
    axs[1].legend(edgecolor="black", facecolor="white", fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, fname))
    plt.close("all")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n_sims",
        type=int,
        default=10000,
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
        default="figures/img/weighting_comps_with_scrambled_scores",
        help="Directory in which generated images should be saved",
    )
    parser.add_argument(
        "--num_scrambled",
        nargs="+",
        default=[0, 500, 900],
        help='List where each element "l" represents an experiment '
        'in which "l" out of 1000 subjects have CATE = 0, on avg',
    )
    args = parser.parse_args()

    for num_scrambled_out_of_1000 in tqdm(args.num_scrambled):
        frac_scrambled = num_scrambled_out_of_1000 / 1000.0
        auc_results = {}
        for method in ["QINI", "AUTOC"]:
            auc_results[method] = []
            for b in range(args.n_sims):
                X, W, Y, Y0, Y1, tau, m_fn = generate_data(
                    N=args.n_samples, p=args.p_treat, seed=b
                )
                scores = get_scores(
                    X, Y, W, e=args.p_treat, m=m_fn, scoring_type="Oracle"
                )
                auc_results[method] += [
                    scaled_RATE(scores[np.argsort(X)], method=method)
                ]

        # Generate a reference set of data with a particular seed
        X, W, _, Y0, Y1, tau, m_fn = generate_data(
            N=args.n_samples, p=args.p_treat, seed=0
        )
        plot_auc_results(
            auc_results,
            X,
            tau,
            Y0,
            Y1,
            save_dir=args.out_dir,
            fname=f"{int(frac_scrambled * 1000)}_out_of_1000_scrambled.png",
        )
