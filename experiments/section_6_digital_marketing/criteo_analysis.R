library(readr)
library(dplyr)

data_loc = "data/criteo-uplift-v2.1.csv"
seed = 20211013
ns = 160000
mtry = 12
honestly = F
min_node_size = 5000

data = read_csv(data_loc)
# print(aggregate(data[,c("visit", "conversion")], list(data$treatment), mean))

set.seed(seed)
subsamp = data %>% mutate(traintest = sample(0:1, n(), replace = TRUE)) %>% group_by(traintest, treatment) %>% sample_n(ns) %>% ungroup() %>% group_by(traintest) %>% group_split()
test_subsamp = subsamp[[2]]
subsamp = subsamp[[1]]
train_features = model.matrix(~.-conversion-treatment-exposure-visit-traintest, data=subsamp)

save(test_subsamp, file="test_subsamp.Rdata")
# load("test_subsamp.Rdata")


print("Average visit:")
print(mean(test_subsamp$visit))
print(t.test(test_subsamp$visit)$"conf.int")
print("Average conversion:")
print(mean(test_subsamp$conversion))
print(t.test(test_subsamp$conversion)$"conf.int")
print("Per arm train")
print(aggregate(subsamp[,c("visit", "conversion")], list(subsamp$treatment), mean))
print("Per arm test")
print(aggregate(test_subsamp[,c("visit", "conversion")], list(test_subsamp$treatment), mean))

test_features = model.matrix(~.-conversion-treatment-exposure-visit-traintest, data=test_subsamp)

model = grf::regression_forest(X=train_features[subsamp$treatment==0, ], Y=subsamp$conversion[subsamp$treatment==0], mtry=mtry, honesty=honestly, min.node.size=min_node_size)
save(model, file="baseline_model.Rdata")
# load("baseline_model.Rdata")
preds = predict(model, newdata=test_features)$predictions

cate_model = grf::causal_forest(X=train_features, Y=subsamp$conversion, W=subsamp$treatment, W.hat=rep(mean(subsamp$treatment), ns*2), mtry=mtry, honesty=honestly, min.node.size=min_node_size)
save(cate_model, file="cate_model.Rdata")
# load("cate_model.Rdata")
preds_cate = predict(cate_model, newdata=test_features)$predictions

y_eval = test_subsamp$visit
eval_model_visit = grf::causal_forest(X=test_features, Y=y_eval, W=test_subsamp$treatment, W.hat=rep(mean(test_subsamp$treatment), ns*2), min.node.size=min_node_size)

y_eval = test_subsamp$conversion
eval_model_conversion = grf::causal_forest(X=test_features, Y=y_eval, W=test_subsamp$treatment, W.hat=rep(mean(test_subsamp$treatment), ns*2), min.node.size=min_node_size)

q = seq(0.001, 1, by = 0.003)
autoc_visit = grf::rank_average_treatment_effect(eval_model_visit, cbind(preds_cate, preds), q=q)
qini_visit = grf::rank_average_treatment_effect(eval_model_visit, cbind(preds_cate, preds), target="QINI", q=q)

autoc_conversion = grf::rank_average_treatment_effect(eval_model_conversion, cbind(preds_cate, preds), q=q)
qini_conversion = grf::rank_average_treatment_effect(eval_model_conversion, cbind(preds_cate, preds), target="QINI", q=q)

save(autoc_visit, autoc_conversion, qini_conversion, qini_visit, file="rates_and_tocs.Rdata")
# load("rates_and_tocs.Rdata")

print("ATE (visit)")
ate = grf::average_treatment_effect(eval_model_visit)
print(ate)
print("ATE (conversion)")
ate = grf::average_treatment_effect(eval_model_conversion)
print(ate)

print("Baseline AUTOC Visit")
print(autoc_visit)
print("Baseline AUTOC Conversion")
print(autoc_conversion)

print("Baseline Qini visit")
print(qini_visit)
print("Baseline Qini Conversion")
print(qini_conversion)
