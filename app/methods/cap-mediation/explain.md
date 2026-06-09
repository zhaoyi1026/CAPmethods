### What CAP-mediation does

**CAP-based mediation analysis** handles a mediation problem where the
**mediator is a subject-level covariance matrix** $\Sigma_i$ (for example, a
brain functional-connectivity matrix), rather than a scalar. A projection
$\theta$ summarizes the mediator as a scalar — its projected log-variance
$\log(\theta^\top \Sigma_i\, \theta)$ — and a hierarchical model links the
exposure, this mediator summary, and the outcome:

$$\text{mediator:}\quad \log(\theta^\top \Sigma_i\, \theta) = \alpha_0 + x_i^\top \alpha + b_i,$$
$$\text{outcome:}\quad Y_i = \gamma_0 + x_i^\top \gamma + \beta\,\log(\theta^\top \Sigma_i\, \theta) + \varepsilon_i,$$

with a subject random intercept $b_i$ in the mediator model. The projection
$\theta$ is **learned from the data**. The quantities of interest are:

- $\alpha$ — exposure → mediator (how the exposure shifts the projected variance),
- $\beta$ — mediator → outcome,
- $\gamma$ — the **direct** effect of the exposure on the outcome,
- **indirect / mediation effect** $\mathrm{IE} = \alpha_{\text{exposure}} \cdot \beta$.

Inference on these effects uses a subject-level **bootstrap**.

### Inputs

- **Mediator** — a list of length $n$; element $i$ is a $T_i \times p$ matrix
  whose covariance is subject $i$'s mediator. Upload as `.rds`/`.RData`, or a long
  CSV (subject-id column + $p$ columns).
- **Exposure / covariates** — an $n \times n_X$ matrix; the **first column is the
  exposure of interest**, any further columns are adjusted covariates. (No
  intercept column — intercepts are estimated internally.)
- **Outcome** — one numeric value per subject (CSV with id + value, or an `.rds`
  numeric vector in subject order).

### Key parameters

- **Number of mediator components** — a fixed count (`nD`) or the DfD criterion.
- **# random initializations** — restarts for the non-convex $\theta$ search.
- **Bootstrap replicates** — subject-level bootstrap for $\alpha$, $\beta$, IE
  inference (0 to skip).

### Outputs

- **Mediation effects** — $\alpha$, $\beta$, direct effect $\gamma$, and the
  indirect effect $\mathrm{IE} = \alpha\cdot\beta$, each with bootstrap SE, test
  statistic, p-value, and confidence interval; plus an effects plot with CIs.
- $\theta$ — the mediator projection loadings.
- A scatter of the outcome vs the mediator score $\log(\theta^\top\Sigma\theta)$,
  colored by exposure (the mediator → outcome path).
- **DfD** across mediator components.

### Built-in example

The example has 100 subjects, a mediator of dimension $p = 10$ measured with
$T_i = 150$ samples each, and a binary treatment. The built-in true indirect
effect is $\mathrm{IE} = \alpha \cdot \beta = 0.8 \times 0.7 = 0.56$. The leading
component recovers the mediator projection $\theta$ (cosine $\approx 1$), the
exposure-on-mediator effect $\hat\alpha \approx 0.8$, and a significant indirect
effect.
