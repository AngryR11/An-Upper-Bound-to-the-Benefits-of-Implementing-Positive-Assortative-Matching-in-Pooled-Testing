# An Upper Bound to the Benefits of Implementing Positive Assortative Matching in Pooled Testing
This project provides a Shiny app that allows users to input relevant parameters—such as disease prevalence, test sensitivity, and specificity, etc.—and instantly obtain the maximum potential benefits that can be achieved by matching together those with similar probabilities of infection when implementing Dorfman Screening.

Details on the methodologies used to compute these upper bounds are available here: [*An Upper Bound to the Benefits of Implementing Positive Assortative Matching in Pooled Testing*](https://ssrn.com/abstract=4779050). Please cite this paper when referencing this app.




## Usage

Download the ``app.R`` and ``function_tight_upper_bound.R`` files and place them on the same directory. Then, open the app.R file and run it on Rstudio and follow these steps:

### Step 1:

Users must first inform the parameters that govern the dilution effect. More precisely, we assume that the probability that a pooled test with $K$ specimens detects an infection, conditional that exactly $I$ of these specimens are infected, is given by the function:
$$
h(I,K)=(1-S_p)+(S_p+S_e-1)\left(\frac{I}{K}\right)^\delta,
$$
where $S_e$ and $S_p$ correspond to the test sensitivity and specificity, respectively, when conducting individual tests. The parameter $\delta\in[0,1]$ governs the dilution effect: the higher higher this parameter, the more the test is subject to dilution effects.

Once $\color{blue}{S_e}$, $\color{blue}{S_p}$ and $\color{blue}{\delta}$ are provided by the user, the app displays the dilution function $h(I,K)$ and proceeds to the next step.

### Step 2:

Afterwards, users are asked to provide the following parameters:


- $\color{blue}{\text{Prevalence}}$: Corresponds to the average probability of infection in the population undergoing testing.

- $\color{blue}{\text{Maximum probability of infection}}$: The maximum probability of infection that can be observed in the population undergoing testing (e.g., for COVID-19, this could be the probability of infection from those who were not vaccinated).

- $\color{blue}{\text{Minimum probability of infection}}$: The lowest probability of infection observed in the population undergoing testing — typically a value close to zero.

- $\color{blue}{\text{Pool size}}$: The number of specimens combined to form a pooled test.

- $\color{blue}{\text{Number of pools}}$: The total number of pools tested per batch. This variable is important because, as shown in [*An Upper Bound to the Benefits of Implementing Positive Assortative Matching in Pooled Testing*](https://ssrn.com/abstract=4779050),  larger batch sizes tend to increase the efficiency gains from ordered pooling. If this entry is left unfilled, the app assumes the number of pools (and therefore the batch size) is arbitrarily large.

- $\color{blue}{\text{Cost of a test (USD)}}$: The cost of performing a single diagnostic test on an individual sample.

## Outputs
The app returns the following as outputs:

- **Maximum reduction in the expected number of tests per subject**: Returns an upper bound to the difference between the expected number of tests obtained when implementing random pooling and the expected number of tests obtained when implementing random pooling, divided by the total number of specimens tested.

- **Maximum cost reduction per subject**:  Corresponds to the previous output multiplied the cost of a test in USD.

- **Maximum reduction in false positives per subject**:  Returns an upper bound to the difference between the expected number of false positives obtained when implementing random pooling and the expected number of false positives obtained when implementing random pooling, divided by the total number of specimens tested.

- **Maximum reduction in false negatives per subject**:  Returns an upper bound to the difference between the expected number of false negatives obtained when implementing random pooling and the expected number of false negatives obtained when implementing random pooling, divided by the total number of specimens tested.

