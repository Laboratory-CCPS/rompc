

Dimensions


| Scalar | |
| --- | --- |
| $n$ | nb. of states |
| $m$ | nb. of controls |
| $N$ | prediction length |
| $r$ | dim. of reduced input space |
| $s$ | $r + 1$ |
| $n_x$ | nb. of inequality conditions on $x$ |
| $n_u$ | nb. of inequality conditions on $u$ |
| $n_f$ | nb. of inequality conditions on $x(T)$ |
| $q$ | total nb. of inequality conditions $q= N n_u + (N - 1) n_x + n_f$ |

| Vector | Dim. | |
| --- | --- | --- |
| $x_k$ | $n$ | states |
| $u_k$ | $m$ | plant inputs |
| $z_k$ | $m$ | inputs w/o contribution of pre-stabilising $K$ |
| $z$ | $Nm$ | $vstack(z_0,\ \dots z_{N-1})$ |
| $z_0$ | $Nm$ | $z_0 = \xi - \tilde z$ |
| $y$ | $q$ | slack variables |
| $\lambda$ | $q$ | lagrange multiplier |
| $\delta$ | $q$ | |
| $p$ | $s$ | |
| $\eta$ | $s$ | |
| $\beta$ | $r + 1 = s$ | |
| $w$ | $N (n + m)$ | control and state sequence $w = vstack(u_0, x_1, \ldots, u_{N-1}, x_N)$ |


| Matrix | Dim. | |
| --- | --- | --- |
| $A$ | $n \times n$ |  |
| $B$ | $n \times m$ |  |
| $Q_{c}$ | $n \times n$ |  |
| $R_{c}$ | $m \times m$ |  |
| $P_{f}$ | $n \times n$ |  |
| $U$ | $Nm \times r$ |  |
| $V$ | $Nm \times (r + 1) = Nm \times s$ | $V = [U\ z_0]$  |
| $W$ | $N (n + m) \times s$ |
| $C$ | $Nn \times N (n + m)$ | |
| $D$ | $Nm \times N (n + m)$ | |
| $E_o$ | $ N (n + m) \times (N (n + m) + s) $ | |
| $Q$ | $(N (n + m) + s) \times (N (n + m) + s)$ | $E_o^T = Q \left[ \begin{array}{c} R \\ 0 \end{array}\right] = [Q_1\ Q_2] \left[ \begin{array}{c} R \\ 0 \end{array}\right] = Q_1 R$ |
| $R$ | $N (n + m) \times N (n + m)$ |
| $Q_1$ | $(N (n + m) + s) \times N (n + m)$ | $Q = [Q_1\ Q_2]$ |
| $Q_2$ | $(N (n + m) + s) \times s$ | |
| $\hat Q_1$ | $(N (n + m) + s) \times Nn$ | $Q = [\hat Q_1\ \hat Q_2\ \hat Q_3]$ |
| $\hat Q_2$ | $(N (n + m) + s) \times Nm + 1$ | |
| $\hat Q_3$ | $(N (n + m) + s) \times r$ | |
| $\check H$ | $(Nm + 1) \times Nm$ | |
| $E$ | $ N (n + m) \times (N (n + m) + r) $ | |
| $\check Q$ | $(N (n + m) + r) \times (N (n + m) + r)$ | $E^T = \check Q \left[ \begin{array}{c} \check R \\ 0 \end{array}\right] = [\check Q_1\ \check Q_2] \left[ \begin{array}{c} \check R \\ 0 \end{array}\right] = \check Q_1 \check R$ |
| $\check Q_1$ | $(N (n + m) + r) \times Nn $ | $\check Q = [\check Q_1\ \check Q_2\ \check Q_3]$ |
| $\check Q_2$ | $(N (n + m) + r) \times Nm $ | |
| $\check Q_3$ | $(N (n + m) + r) \times r $ | |
| $\check R$ | $N (n + m) \times N (n + m)$ |  |
| $\check R_{11}$ | $Nn \times Nn$ | $\check R = \left[ \begin{array}{cc} \check R_{11} & \check R_{12} \\ 0 & \check R_{22} \end{array} \right]$ |
| $\check R_{12}$ | $Nn \times Nm$ |  |
| $\check R_{22}$ | $Nm \times Nm$ |  |
| $\check R$ | $N (n + m) \times N (n + m)$ |  |
| $\tilde R$ | $Nm \times Nm$ | |