# Lab 13
## Wakefield exitation

In the recent lecture, we established the nonlinear wave equation for the 
vector potential A of a weakly relativistic, linearly polarized 
laser pulse in a plasma of (normalized) constant density *n*
is

<p align="center">
<img src="stuffy_stuff/eq2.png" width="300">
</p>



However, the above equation assumes that the plasma density is constant, which
is not correct for high intensity laser pulses. Via its ponderomotive force the
laser pulse acts on the (electron) plasma density. The coupling between light wave and plasma density distibution can be modeled by the following system,
where <img src="stuffy_stuff/eq3.png" width="75">.

<p align="center">
<img src="stuffy_stuff/eq1.png" width="340">
</p>

This set of coupled PDEs can be solved via a splitting method. The idea is to propagate *A* according to the first equation for one time-step (and keeping *dn*
constant) and then propagate *dn* for one time step using (and keeping *A* constant). Once the cylce is complete, we advanced the full system *(A,dn)* for one time-step.

Second-order finite differences lead to the explicit equations

<p align="center">
<img src="stuffy_stuff/eq4.png" width="500">
</p>
<p align="center">
<img src="stuffy_stuff/eq5.png" width="500">
</p>

---
### Completing the code

Implement two function, `stepA` and `stepdn` that use the above equations
to determine <img src="stuffy_stuff/eq6.png" width="45">
and <img src="stuffy_stuff/eq7.png" width="45">

---
### Initial conditions

For the initial conditions, we use a Gaussian laser pulse
<p align="center">
<img src="stuffy_stuff/eq8.png" width="400">
</p>
and assume <img src="stuffy_stuff/eq9.png" width="50">. Even though these initial
conditions are not consistent with the set of equations they will be good 
enough for our purposes.