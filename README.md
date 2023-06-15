# Smooth IDR
Extension of Isotonic Distributional Regression (https://github.com/AlexanderHenzi/isodistrreg) to its kernel-smoothed version, which yields a smooth and consistent estimator of the conditional CDF under the assumption of stochastic monotonicity. This also provides a consistent estimator of the conditional density, with two bw selection approaches.

## Introduction

Isotonic distributional regression (IDR) aims at modeling the conditional CDF of a random variable $Y$ given $X$ under the assumption that $Y$ is in a monotonic relationship with $X.$ Formally, we assume the conditional CDF $F_x(y)$ defined for $(x,y) \in \mathcal{X}\times \mathbb{R}$ to be non-increasing in $x$ for any $y \in \mathbb{R}.$
The IDR estimator is obtained by minimizing the Continuous Ranked Probability Score (CRPS).
When the observations of $X$ are given by $x_1 < \dots < x_n$ and the corresponding observations of $Y$ by $y_1,\dots,y_n,$ the IDR conditional CDF of $Y$ at a threshold $y \in \mathbb{R}$ given $x_i$ is given by

$$\begin{equation}
    H(x_i,y)=\min_{\substack{j:j\leq i}} \max_{\substack{h:h\geq i}} \frac{1}{h-j+1}\sum_{k=j}^h \mathbb{1}_{[Y_k\leq y]}.
\end{equation}$$

For values of $x \in \mathcal{X},$ the IDR CDF at a threshold $y \in \mathbb{R}$ is defined as a linear interpolation between $H(x_j,y)$ and $H(x_{j+1},y),$ where $1\leq j\leq n$ is such that 
$x_j < x< x_{j+1}.$ 
The function $H(x,\cdot)$ is a piecewise continuous function with jumps at the unique observed values of $Y$ given by $\widetilde{y}_1,\dots,\widetilde{y}_m.$


We obtain a smooth estimator of the conditional CDF using kernel-smoothing. We define the kernel-smoothed IDR CDF as 

$$\begin{equation}
        \widetilde{H}(x,y)=\int_\mathbb{R}H(x,t)K_h(y-t) \ dt,
\end{equation}$$
    
where $K_h(x)=\kappa(x/h)/h$ and $\kappa$ is a valid kernel function and $h>0$ is the bandwidth.

Since $H(x,y)$ is discrete with jumps at $\widetilde{y}_i$, we can write

```math
\begin{equation}
\widetilde{H}(x,y) = \sum_{j=1}^m H(x,\widetilde{y}_j) \int_{\widetilde{y}_j}^{\widetilde{y}_{j+1}} K_h(y-t) \ dt.
\end{equation}
````

where we define the boundary cases $\widetilde{y}_{m+1} = \infty$ and $\widetilde{y}_0 = -\infty.$

From the last equation, we can infer on the condtional density of $Y|X=x$ by computing the derivative of the smooth IDR:
```math
$$\begin{align} \widetilde{h}(x,y)&= \frac{\partial}{\partial y}\widetilde{H}(x,y)= \sum_{j=1}^mH(x,\widetilde{y}_j)(K_h(y-\widetilde{y}_j)-K_h(y-\widetilde{y}_{j+1})) \notag \\
    &=\sum_{j=1}^mH(x,\widetilde{y}_j)K_h(y-\widetilde{y}_j)-\sum_{j=1}^mH(x,\widetilde{y}_j)K_h(y-\widetilde{y}_{j+1})\notag \\
    &= \sum_{j=1}^mw_j(x)K_h(y-\widetilde{y}_j),
\end{align}$$
````
where
```math 
w_j(x) = H(x,\widetilde{y}_j) - H(\widetilde{y}_{j-1}).
````

It can be shown under given assumptions, that the IDR conditional CDF and its smooth version are consistent estimators of the conditional CDF. Furthermore under additional regularity assumptions, also the conditional IDR density is a consistent estimator of the conditional density.

## Performance assessment and bandwidth selection

The smooth IDR CDF and density require tuning of the bandwith $h,$ for which we provide two possible approaches.
One possible approach is to compute the cross-validation error for different values of $h$ and choose the value that minimizes the score. Using the $\mathrm{logS}$ scoring rule, the leave-one-out cross-validation error is defined by 

```math
\begin{equation*}
    \mathrm{CV}(h):=-\frac{1}{n}\sum_{i=1}^n \log{\widetilde{h}_{x_{-i}}(x_i,y_i)}.
\end{equation*}
```` 

The IDR CDF $H$ can be efficiently computed with the PAV algorithm in $O(log(n)n)$ time. However, computing $\mathrm{CV}(h)$
requires $O(mlog(n)n^2)$ computational steps. The cross-validation error can be estimated using a new procedure called 'one-fit grid search' which requires only one instead of n model generations. 
The criterion is defined as

$$\begin{equation*}
    \mathrm{OF}(h):=-\frac{1}{n}\sum_{i=1}^n \log{\bar{h}_{-i}(x_i,y_i)},
\end{equation*}$$


where $\log{\bar{h}_{-i}(x_i,y_i)}$ is given by setting the $i^{th}$ weight to zero and re-scaling the other weights of the IDR density. 
Computing $\mathrm{OF}(h)$ reduces the computational effort to $O(nmlog(n))$ steps.

The second approach is a local bandwidth estimator which is obtained by minimizing an upper bound of the $\mathrm{MSE}$ which converges to zero as $n$ approaches infinity. The estimator is defined as:

 ```math
 \begin{equation}
    \hat{h}^*_n(y)=c\bigg(\frac{4\kappa(0)\big(\frac{\log(n)}{n}\big)^{1/3}}{\big|\widetilde{h}''_x(y)\big|\int_\mathbb{R} u^2\kappa(u)\ du}\bigg)^\frac{1}{3}.
\end{equation}
````


The constant $c$ is unknown in general and has to be chosen a priori or by minimizing $\mathrm{OF}$ or $\mathrm{CV}.$

It is thus advantageous, that the absolute difference $\max_{c \in [c_{\min}, c_{\max}]}|\mathrm{OF}(h_{c,n}) - \mathrm{CV}(h_{c,n})|$ converges in probability to zero, for $h_{c,n}=(\log(n)/n)^{1/9}c$ and  $0< c_{\min} \leq c_{\max} < \infty.$

## Experiments

To visualize the advantage of $\mathrm{OF}(h),$ we simulate for different sample sizes the random variable 
```math
Y|X\sim \textup{Gamma}(\textup{shape}=X,\textup{scale}=1/\sqrt{X}).
````
To furthermore illustrate the consistent estimation of $\mathrm{CV}(h),$ we choose for every $n$ the bandwidth $h_n= (\log{(n)}/n)^{1/9}.$ When comparing $\mathrm{CV}(h_n)$ and $\mathrm{OF}(h_n),$ we observe that $\mathrm{OF}$ closely matches $\mathrm{CV}$ while being much faster to compute.





  
