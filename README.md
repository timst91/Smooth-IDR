# Smooth IDR
Extension of Isotonic Distributional Regression (https://github.com/AlexanderHenzi/isodistrreg) to its kernel-smoothed version, which yields a smooth and consistent estimator of the conditional CDF under the assumption of stochastic monotonicity. This also provides a consistent estimator of the conditional density, with two bw selection approaches.

## Introduction

Isotonic distributional regression (IDR) aims at modeling the conditional CDF of a random variable $Y$ given $X$ under the assumption that $Y$ is in a monotonic relationship with $X.$ Formally speaking, this is equivalent to the conditional CDF $F_x(y)$ defined for $(x,y) \in \mathcal{X}\times \mathbb{R}$ is non-increasing in $x$ for any $y \in \mathbb{R}.$
When the observations of $X$ are given by $x_1<\dots<x_n$ and the corresponding observations of $Y$ by $y_1,\dots,y_n,$ the IDR conditional CDF of $Y$ at a threshold $y \in \mathbb{R}$ given $x_i$ is defined as 

$\begin{equation}\label{eq:IDR}
    H(X_i,y)=\min_{\substack{j:j\leq i}} \max_{\substack{h:h\geq i}} \frac{1}{h-j+1}\sum_{k=j}^h \mathds{1}_{[Y_k\leq y]}.
\end{equation}$

For values of $x \in \mathcal{X},$ the IDR CDF at a threshold $y \in \mathbb{R}$ is defined as a linear interpolation between $H(x_j,y)$ and $H(x_{j+1},y),$ where $1\leq j\leq n$ is such that 
$x_j <x<x_{j+1}.$ 
The function $H(x,\cdot)$ is a piecewise continuous function with jumps at the unique observed values of $Y$ given by $\Tilde{y}_1,\dots,\Tilde{y}_m.$

We obtain a smooth estimator of the conditional CDF using kernel-smoothing. We define the kernel-smoothed IDR CDF as 
$\begin{equation}
        \Tilde{H}(x,y)=\int_\mathbb{R}H(x,t)K_h(y-t) \ dt,
    \end{equation}$
where $K_h(x)=\kappa(x/h)/h$ and $\kappa$ is a valid kernel function and $h>0$ is the bandwidth.

Since $H(x,y)$ is discrete with jumps at $\Tilde{y}_i,$ we can write $\begin{equation}
    \Tilde{H}(x,y)=\sum_{j=1}^mH(x,\Tilde{y}_j)\int_{\Tilde{y}_j}^{\Tilde{y}_{j+1}}K_h(y-t) \ dt,
\end{equation}$
where we define the boundary cases $\tilde{y}_{m+1} = \infty$ and $\tilde{y}_0 = -\infty.$

From the last equation, we can infer on the condtional density of $Y|X=x$ by computing the derivative of the smooth IDR:
$\begin{align}\label{eq: weighted kernell form of IDR density}
    \Tilde{h}(x,y)&= \frac{\partial}{\partial y}\Tilde{H}(x,y)= \sum_{j=1}^mH(x,\Tilde{y}_j)(K_h(y-\Tilde{y}_j)-K_h(y-\Tilde{y}_{j+1})) \notag \\
    &=\sum_{j=1}^mH(x,\Tilde{y}_j)K_h(y-\Tilde{y}_j)-\sum_{j=1}^mH(x,\Tilde{y}_j)K_h(y-\Tilde{y}_{j+1})\notag \\
    &= \sum_{j=1}^mw_j(x)K_h(y-\Tilde{y}_j),
\end{align}$
where $w_j(x)=H(x,\Tilde{y}_j)-H(x,\Tilde{y}_{j-1}).$

It can be shown under given assumptions, that the IDR conditional CDF and its smooth version are consistent estimators of the conditional CDF. Furthermore under additional regularity assumptions, also the conditional IDR density is a consistent estimator of the conditional density.

##



  