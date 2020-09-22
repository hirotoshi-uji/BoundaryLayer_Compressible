# 圧縮性境界層方程式の数値解

圧縮性流体における境界層方程式を解いた．

境界層方程式を変数変換により，常微分方程式に直し，shooting method により解いた．4次精度のRunge-Kutta法を用いた．  
- $f^{\prime \prime \prime} + \frac{1}{2} ff^{\prime \prime} = 0$  
- $g^{\prime \prime} + \frac{Pr}{2}g^{\prime}f+Pr(\gamma -1)Ma^2f^{\prime \prime}{}^2 = 0$  
- $f = f^{\prime} = 0$ at $\eta = 0$  
- $f^{\prime} \rightarrow 1$ at $\eta \rightarrow \infty$
- $g^{\prime} = 0 $ at $\eta=0$  
- $g \rightarrow 1$ at $\eta \rightarrow 1$  

