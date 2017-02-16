foo <-
structure(list(structure(list(structure("findpars", Rd_tag = "VERB")), Rd_tag = "\\name"), 
    structure("\n", Rd_tag = "TEXT"), structure(list(structure("1.1", Rd_tag = "VERB")), Rd_tag = "\\Rdversion"), 
    structure("\n", Rd_tag = "TEXT"), structure(list(structure("findpars", Rd_tag = "VERB")), Rd_tag = "\\alias"), 
    structure("\n", Rd_tag = "TEXT"), structure(list(structure("\n", Rd_tag = "TEXT"), 
        structure("Fit and compare mortality models.\n", Rd_tag = "TEXT")), Rd_tag = "\\title"), 
    structure("\n", Rd_tag = "TEXT"), structure(list(structure("\n", Rd_tag = "TEXT"), 
        structure("Fit and compare all supported models to a dataset evaluate the likelihood of a parameter differing between two datasets. Instead of using 'reasonable' starting values, this function takes advantage of the fact that many mortality models are nested in each other and all are ultimately nested in the exponential model whose sole parameter can be calculated directly. Each successively more complex model uses the parameters fitted to the model it is nested in as starting values.\n", Rd_tag = "TEXT")), Rd_tag = "\\description"), 
    structure("\n", Rd_tag = "TEXT"), structure(list(structure("\n", Rd_tag = "RCODE"), 
        structure("findpars(x, y = NULL, cx = NULL, cy = NULL, nil = 0, bnil = 0, wbnil = 1, pf = mean, label = NULL, summary = F, id = 0, tlog = F, digits = 22, sig = 0.05, models = NULL)\n", Rd_tag = "RCODE")), Rd_tag = "\\usage"), 
    structure("\n", Rd_tag = "TEXT"), structure(list(structure("\n", Rd_tag = "TEXT"), 
        structure("  ", Rd_tag = "TEXT"), structure(list(list(
            structure("x", Rd_tag = "TEXT")), list(structure("\n", Rd_tag = "TEXT"), 
            structure("A vector of ages at death (integers, usually representing days).\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("y", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("An optional second vector of ages at death (integers, usually representing days).\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("cx", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("A vector of ones and zeros, with ones representing a natural death event, and zeros representing a censored event. If omitted, it is assumed that all the deaths are natural in the control group.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("cy", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("A vector of ones and zeros, with ones representing a natural death event, and zeros representing a censored event. If omitted, it is assumed that all the deaths are natural in the experimental group.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("nil", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("A starting numeric value for a newly added parameter. For example, the Logistic model uses the best-fit parameters for the Gompertz model as its starting values, but the Gompertz model has two parameters while the Logistic model has three. The third parameter takes on the value of the nil argument.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("bnil", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("Like nil, but bnil only gets substituted into the b parameter.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("wbnil", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("The numeric value to substitute for missing values of the 'b' parameter when finding parameters for the Weibull model. Defaults to 1.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("pf", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("A function to call when two different parameters are constrained and need to produce a single starting value. In addition to ", Rd_tag = "TEXT"), 
                structure(list(structure("mean", Rd_tag = "RCODE")), Rd_tag = "\\code"), 
                structure(" (default), other valid functions include ", Rd_tag = "TEXT"), 
                structure(list(structure("median", Rd_tag = "RCODE")), Rd_tag = "\\code"), 
                structure(", ", Rd_tag = "TEXT"), structure(list(
                  structure("gmean", Rd_tag = "RCODE")), Rd_tag = "\\code"), 
                structure(", ", Rd_tag = "TEXT"), structure(list(
                  structure("max", Rd_tag = "RCODE")), Rd_tag = "\\code"), 
                structure(", and ", Rd_tag = "TEXT"), structure(list(
                  structure("min", Rd_tag = "RCODE")), Rd_tag = "\\code"), 
                structure(".\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("label", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("A character prefix based on which the output files are named. You may specify a file-path in the prefix, but do not specify a file extension as the extension '.txt' will automatically be appended to the filenames. If the label is not specified, findpars will not save an output file.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("summary", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("Logical value telling findpars whether or not to echo a summary table to the console.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("id", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("A numeric value, copied to the 'id' column of the output table. This is used to distinguish different iterations when findpars is called by another function or script.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("tlog", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("Logical value passed as the tlog argument to the opsurv function, instructing it to find best fits for the natural log of the starting values and then convert them back before outputting them. Maybe result in improved speed but decreased accuracy and is still being tested.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("digits", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("Integer specifying many digits to show in the summary output to the console. Only useful if summary is set to TRUE.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("sig", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("Significance cutoff used by the chi-squared test.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure("  ", Rd_tag = "TEXT"), 
        structure(list(list(structure("models", Rd_tag = "TEXT")), 
            list(structure("\n", Rd_tag = "TEXT"), structure("Optional character vector containing the names of models to test. If you specify only the models you're interested in, findpars is smart enough to automatically fill in the prerequisite models and fit them as well.\n", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT")), Rd_tag = "\\arguments"), 
    structure("\n", Rd_tag = "TEXT"), structure(list(structure("\n", Rd_tag = "TEXT"), 
        structure("This function takes one or two vectors of ages at death, fits various mortality models to them, and compares the fits of each mortality model to the next simpler mortality model that it's based on to find the best models. If two vectors (x and y) have been given, for each parameter of each joint model, this function will compare the unconstrained joint model to a joint model that is constrained such that that parameter is the same between the two groups. A significantly worse fit than the unconstrained model indicates the parameter being tested is in fact different. While running the function draws a progress bar on the console and prints the name of a model (the unconstrained version of it) each time a fit is complete. This is normal behavior whose purpose is to reassure the user that the function hasn't hung and is indeed working, because for some datasets it may take a long time to converge on the best parameters.\n", Rd_tag = "TEXT")), Rd_tag = "\\details"), 
    structure("\n", Rd_tag = "TEXT"), structure(list(structure("\n", Rd_tag = "TEXT"), 
        structure("A data.frame, invisibly returned. If summary=T, the table is also printed to the console. If a label argument is given, the data is also saved to a tab-delimited file of that name. Cells in the table for which there are no applicable values are filled with NAs. The output contains the following columns:\n", Rd_tag = "TEXT"), 
        structure(list(list(structure("MLE", Rd_tag = "TEXT")), 
            list(structure("Maximum likelihood estimate for the model described in that row.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("# cycles", Rd_tag = "TEXT")), list(structure("The number of iterations the optimization function had to go through before it converged.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("OK?", Rd_tag = "TEXT")), list(structure("If this model has a better fit than the next simplest model, this value will be '1'. If the fit is the same, this value will be 0. If the fit is worse, this value will be -1 and in such cases you are encouraged to send your data and your command history to the author because you have most likely uncovered a bug in the software.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("a1", Rd_tag = "TEXT")), list(structure("The lambda parameter of the first group, for all models.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("b1", Rd_tag = "TEXT")), list(structure("The beta parameter of the first group for a Weibull model and the gamma parameter of the first group for a Gompertz, Logistic, Gompertz-Makeham, or Logistic-Makeham model.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("c1", Rd_tag = "TEXT")), list(structure("The c parameter of the first group for a Gompertz-Makeham or Logistic-Makeham model.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("s1", Rd_tag = "TEXT")), list(structure("The s parameter of the first group for a Logistic or Logistic-Makeham model.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("a1", Rd_tag = "TEXT")), list(structure("The lambda parameter of the second group, for all models.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("b1", Rd_tag = "TEXT")), list(structure("The beta parameter of the second group for a Weibull model and the gamma parameter of the second group for a Gompertz, Logistic, Gompertz-Makeham, or Logistic-Makeham model.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("c1", Rd_tag = "TEXT")), list(structure("The c parameter of the second group for a Gompertz-Makeham or Logistic-Makeham model.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("s1", Rd_tag = "TEXT")), list(structure("The s parameter of the second group for a Logistic or Logistic-Makeham model.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("model", Rd_tag = "TEXT")), list(structure("The model being fitted in this row.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("id", Rd_tag = "TEXT")), list(structure("A copy of the 'id' argument.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("nil", Rd_tag = "TEXT")), list(structure("A copy of the 'nil' agrument.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("bnil", Rd_tag = "TEXT")), list(structure("A copy of the 'bnil' argument.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("LR", Rd_tag = "TEXT")), list(structure("Log ratio for this fit-- i.e. log(MLE this model) - log(MLE next simpler model).", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("AIC", Rd_tag = "TEXT")), list(structure("Akaike's Information Criterion for this model fit. This can be used for comparing the fit of models that are not nested in each other.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("BIC", Rd_tag = "TEXT")), list(structure("Bayesian Information Criterion for this model fit, which is similar to the AIC but with the additional feature of penalizing for the number of parameters.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("p (chi squared)", Rd_tag = "TEXT")), list(
            structure("The significance level of the observed LR assuming a chi-squared null distribution of LRs with one degree of freedom.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT"), structure(list(list(
            structure("sig?", Rd_tag = "TEXT")), list(structure("Whether or not the significance level is below the threshold set by the 'sig' argument.", Rd_tag = "TEXT"))), Rd_tag = "\\item"), 
        structure("\n", Rd_tag = "TEXT")), Rd_tag = "\\value"), 
    structure("\n", Rd_tag = "TEXT"), structure(list(structure("\n", Rd_tag = "TEXT"), 
        structure("Pletcher,S.D., Khazaeli,A.A., and Curtsinger,J.W. (2000). Why do life spans differ? Partitioning mean longevity differences in terms of age-specific mortality parameters. Journals of Gerontology Series A-Biological Sciences and Medical Sciences 55, B381-B389\n", Rd_tag = "TEXT")), Rd_tag = "\\references"), 
    structure("\n", Rd_tag = "TEXT"), structure(list(structure("\n", Rd_tag = "TEXT"), 
        structure("Alex F. Bokov\n", Rd_tag = "TEXT")), Rd_tag = "\\author"), 
    structure("\n", Rd_tag = "TEXT"), structure("\n", Rd_tag = "TEXT"), 
    structure(list(structure("\n", Rd_tag = "TEXT"), structure(list(
        structure(list(structure("opsurv", Rd_tag = "TEXT")), Rd_tag = "\\link")), Rd_tag = "\\code"), 
        structure(", ", Rd_tag = "TEXT"), structure(list(structure(list(
            structure("simdist", Rd_tag = "TEXT")), Rd_tag = "\\link")), Rd_tag = "\\code"), 
        structure(", ", Rd_tag = "TEXT"), structure(list(structure(list(
            structure("empdist", Rd_tag = "TEXT")), Rd_tag = "\\link")), Rd_tag = "\\code"), 
        structure("\n", Rd_tag = "TEXT")), Rd_tag = "\\seealso"), 
    structure("\n", Rd_tag = "TEXT"), structure("\n", Rd_tag = "TEXT"), 
    structure(list(structure(" survival ", Rd_tag = "TEXT")), Rd_tag = "\\keyword"), 
    structure("\n", Rd_tag = "TEXT")), class = "Rd")
