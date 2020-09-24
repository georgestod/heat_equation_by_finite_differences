# Linear heat equation solver by finite differences

The heat equation is a beautiful equation that can be derived from first principles. Stating the energy conservation law leads to,
\begin{equation*}
\frac{d}{dt} \Big( \underbrace{\int_D \rho \cdot c_p \cdot T \cdot dv}_{\text{total energy}} \Big ) =  \underbrace{\int_{\partial D} k \nabla \cdot T \cdot n \cdot dS}_{\text{power through boundaries}} + \underbrace{\int_D q \cdot dv}_{\text{power source term}}
\end{equation*}
which becomes by introducing the divergence operator,
\begin{equation*}
\frac{d}{dt} \Big( \int_D \rho \cdot c_p \cdot T \cdot dv \Big ) =  \int_{D} \nabla \cdot (k \nabla \cdot T) dv + \int_D q \cdot dv
\end{equation*}
which means for every point of the domain $D$,
\begin{equation*}
\frac{\partial (\rho \cdot c_p \cdot T)}{\partial t} = \nabla \cdot (k \cdot \nabla T) + q
\end{equation*}
If we assume the density and specific heat of the $D$ are time-independent,
\begin{equation*}
\rho \cdot c_p \cdot \frac{\partial T}{\partial t} = \nabla \cdot (k \cdot \nabla T) + q
\end{equation*}
Now, if we add the assumption that isotropy (same conductivity in all directions) on $D$,
\begin{equation*}
\rho \cdot c_p \cdot \frac{\partial T}{\partial t} = k\cdot \nabla^2 T + q
\end{equation*}
Finally, if we limit ourselves to 2 dimensions and suppose there is no energy input in our system,
\begin{equation}
\frac{\partial T}{\partial t} = \frac{k}{\rho \cdot c_p} \Big( \frac{\partial^2 T}{\partial x^2} +\frac{\partial^2 T}{\partial y^2} \Big)
\label{eq:THE_EQUATION}
\end{equation}
\subsection{Discretizing time and space}
In its simplest form, finite differences consist in discretizing the \textit{time} and \textit{space} domains as follows,
\begin{equation}
\Big(\ \frac{\partial T}{\partial t} \Big)_{i,j}^{n} \approx \frac{T_{i,j}^{n+1}-T_{i,j}^{n}}{\Delta t}
\end{equation}
%
\begin{equation}
\Big(\ \frac{\partial^2 T}{\partial x^2} \Big)_{i,j}^{n} \approx \frac{T_{i+1,j}^{n}-2 T_{i,j}^{n}+T_{i-1,j}^{n} }{\Delta x^2}
\end{equation}
%
\begin{equation}
\Big(\ \frac{\partial^2 T}{\partial y^2} \Big)_{i,j}^{n} \approx \frac{T_{i,j+1}^{n}-2 T_{i,j}^{n}+T_{i,j-1}^{n} }{\Delta y^2}
\end{equation}
%
where $n$ is a component of the discretized time vector, $(i,j)$ are the coordinates of a point in the $xy$ plane, $\Delta t$, the time step and $\Delta x$ and $\Delta y$ the space between points on the desired computation domain.
By combining the expressions above and the equation \ref{eq:THE_EQUATION}, the temperature can be computed on every point of the computation domain as a function of the previous time step,
\begin{equation}
T_{i,j}^{n+1} = T_{i,j}^{n}+ \Delta t \frac{k}{\rho \cdot c_p} \Big(\frac{T_{i+1,j}^{n}-2 T_{i,j}^{n}+T_{i-1,j}^{n} }{\Delta x^2} + \frac{T_{i,j+1}^{n}-2 T_{i,j}^{n}+T_{i,j-1}^{n} }{\Delta y^2}\Big)
\end{equation}
\subsection{Numerically stable solutions}
The previous equation can be recasted in the following matrix form,
\begin{equation}
T^{n+1} = \Big( \frac{\alpha \Delta t}{h^2}\cdot A+I_n\Big)\cdot T_n + \frac{\alpha \Delta t}{h^2}\cdot B
\label{eq:suite}
\end{equation}
where $B$ encapsulates static boundary conditions (as the ones used in this work), $h=\Delta x=\Delta y$ and $\alpha = \frac{k}{\rho \cdot c_p}$ is the \textit{diffusion} coefficient.  To reduce the computational cost, it is recommended when using this formulation to specify in MATLAB or equivalent that $A$ is a sparse matrix.

The expression \ref{eq:suite} is a finite arithmetic progression, so to make sure it will converge, nonnegativity is required,
\begin{equation}
 \frac{\alpha \Delta t}{h^2}\cdot A+I_n \ge 0
\end{equation}
By looking at the smallest terms, which are in the diagonal, we get the largest time step allowed for numerical stability,
\begin{equation}
\Delta t \leq \frac{h^2}{\gamma \cdot \alpha}
\end{equation}
where $\gamma = 4$. All the previous steps can be generalized to the 3D case, one would find $\gamma = 6$.
