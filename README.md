# OL-CCS
Open-Loop Coherent Compressed Sensor

Author: M. D. Sudeera H. Gunathilaka

Github: https://github.com/SuhiG

# Installation

## local install (windows)

1) create a python environment (e.g. CCS)
2) clone the python repository found at https://github.com/NTTRI-PHI-Algorithms/CCS
3) install via "pip install ." at the root of the repository
4) run "python ./CSS/random.py",
"python ./CSS/random_observation_noise.py"


# 1) Description

## Introduction

Compressed sensing (CS) is a method of reconstructing a high-dimensional signal or image based on highly downsampled measurements.

The Coherent Compressed Sensor (CCS) is designed to solve the problem of $l_0$-Regularised Compressed Sensing, which is a combinatorial optimisation problem. 

## Model

The problem Hamiltonian can be stated as follows.

$$ H = \sum_{r<r\'}^{N} \sum_{k = 1}^{M} $$

The above equation shows an observed signal $y \in \mathbb{R}^M$ and an observation matrix $A \in \mathbb{R}^{M\times N}$. The term $R \in \mathbb{R}^N$ and $\sigma \in \left\{{0,1}\right\}^N$ correspond to the source signal and support vector, respectively.
In this case, the quadratic optimisation part of the problem (signal - $R_{r}$) is done by using the Classical Digital Processor (CDP) while the combinatorial optimisation part of the problem (support - $\sigma_{r}$) is performed by the CIM in alternate steps.

For the CIM, the injection and local fields can be specified as follows.

$$    \left(\dfrac{dc_{r}}{dt}\right)_{inj,r} = \left(|h_r| - \eta\right). $$

$$h_{r} = -{\sum_{r' = 1 (\neq r)}^{N}\sum_{k = 1}^{M}} A_r^k A_{r'}^k R_{r'}H(c_{r'}) + \sum_{k=1}^M A_{r}^k y^{k} . $$

$H(c_r)$ is the binarised in-phase amplitude by the Heaviside step function. $\eta$ is the threshold which is related to the $l_0$-regularisation parameter $\lambda$ by $\eta = \sqrt{2\lambda}$ according to the Maxwell rule.
Open-Loop CCS is composed of the following stochastic differential equations.

$$
        \frac{d}{dt}c_r = \left[-1 + p - {\left(c_r^2 + s_r^2\right)} \right]c_r + \widetilde{K}\left(\dfrac{dc_{r}}{dt}\right)_{inj,r} + {g^2}\sqrt{\left(c_r^2 + s_r^2\right) + \frac{1}{2}} W_{r,1},$$


$$\frac{d}{dt}s_r = \left[-1 - p - {\left(c_r^2 + s_r^2\right)}\right]s_r + {g^2}\sqrt{\left(c_r^2 + s_r^2\right) + \frac{1}{2}} W_{r,2} .
$$

Here, in-phase and quadrature-phase normalised amplitudes are represented as $c_r$ and $s_r$ respectively. $p$ is the normalised pump rate. $W_{r,1}$ and $W_{r,2}$ are independent real Gaussian noise processes. The term ${g^2}$ indicates the saturation parameter while $\widetilde{K}$ indicates the normalised feedback strength. 

In the CDP, we solve the following system of equations using the Jacobi method.

$$ R_{r}\sum_{k = 1}^{M} \left(A_{r}^{k}\right)^2 = \sigma_{r}\mathbb{H}_{r}, $$


$$ \mathbb{H}_{r} = -\sum_{r' = 1 (\neq r)}^{N}\sum_{k = 1}^{M} A_{r}^{k}A_{r'}^{k}R_{r'}\sigma_{r'} + \sum_{k =1}^{M} A_{r}^{k}y^{k} .$$


## Pseudo-code

    FOR rep IN 0:reps
        FOR t IN 0:tmax
            Set h_i to J*(rCDP*sigCIM) + hz                      #Matrix-Vector multiplication
            Set c to c + (-(1 - p + c**2 + s**2)*c + JJ * (abs(h_i) - eta) + g * sqrt(c**2 + s**2 + 0.5) * w1)*dt
            Set s to s + (-(1 + p + c**2 + s**2)*s + g * sqrt(c**2 + s**2 + 0.5) * w2)*dt
            Set sigCIM to ((c/abs(c))+1)/2
        END FOR
        FOR tt IN 0:ttmax
            Set hCDP to J*(rCDP*sigCIM) + hz                     #Matrix-Vector multiplication
            Set rCDP to rCDP + (-rCDP + hCDP)*dt1
        END FOR
    END FOR
	
## Benchmark

OL-CCS performance is evaluated for randomly generated CS problem instances and Magnetic Resonance Imaging (MRI) of a brain. The performance accuracy is evaluated by calculating the Root-Mean-Squared-Error (RMSE) as $$ \sqrt{\dfrac{1}{N} \sum_{r=1}^N \left(R_r\sigma_r - x_r\xi_r\right)^2}, $$

where $R_r$ and $\sigma_r$ is the estimated signal and support from OL-CCS and $x_r$ and $\xi_r$ is the correct signal and support for the CS problem instance. 

## References

[1]. T. Aonishi, K. Mimura, M. Okada, and Y. Yamamoto, “L0 regularization based compressed sensing with quantum–classical hybrid approach,” Quantum Science and Technology 7, 035013 (2022)

[2]. M. D. S. H. Gunathilaka, S. Kako, Y. Inui, K. Mimura, M. Okada, Y. Yamamoto, and T. Aonishi, “Effective implementation of L0-regularised compressed sensing with chaotic-amplitude-controlled coherent ising machines,” Scientific Reports 13, 16140 (2023).
