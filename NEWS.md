# brm 1.1
    - Fixed a bug when there is only an intercept term for alpha or beta
	- Introduced vectorised functions for getProbRR and getProbRD which makes the model fitting 10-100 times faster
	- Fixed a bug under Linux: matrix objects now also inherit from class "array", namely, e.g., class(diag(1)) is c("matrix", "array") which invalidates code assuming that length(class(obj)) == 1.
	- Use the "one-step" doubly robust estimator as in Richardson, Robins and Wang (2017) eqn. (3.4) - runs up to 50 times faster
