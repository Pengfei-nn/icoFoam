## icoFoam的IO机制介绍

```C++
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
runTime.loop()
runTime.write()
```

"setRootCaseLists.H"：#include "listOptions.H"  #include "setRootCase.H"  #include "listOutput.H"

"createTime.H"：定义时间类  runTime   Foam::Time runTime(Foam::Time::controlDictName, args);

...

- 对象注册机的介绍

  ...

## icoFoam中头文件介绍

```c++
#include "fvCFD.H"
#include "pisoControl.H"

#include "postProcess.H

#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"

#include "createFields.H"
#include "initContinuityErrs.H"

 #include "continuityErrs.H"
```

"fvCFD.H" ：有限体积法头文件，其中包含了 fvMesh.H，fvMatrix.H等。

"pisoControl.H" ：pisoControl类的声明

...



## PISO 算法流程以及在icoFoam中的实现

$$
\left\{\begin{array}{l}\nabla \cdot\boldsymbol{U}=0 \\\frac{\partial \boldsymbol{U}}{\partial t}+\nabla\cdot(\boldsymbol{U}\otimes \boldsymbol{U}) - \nabla\cdot(\nu\nabla \boldsymbol{U})= -\nabla P\end{array}\right.
$$

- 问题的求解难点在于动量的非线性项，以及动量，压力的耦合。

### 动量预测

对动量方程在控制体上积分，并利用高斯散度定理得到：
$$
\int_{V}\frac{\partial \boldsymbol{U}}{\partial t}dV + \oint_{S}\boldsymbol{U} \boldsymbol{U}\cdot d\boldsymbol{S} -  \oint_{S}\nu\nabla\boldsymbol{U}\cdot d\boldsymbol{S} = \oint_{S}-Pd\boldsymbol{S}
$$
记面$f$上的通量$\Phi_f = \boldsymbol{U}_f \cdot d\boldsymbol{S}$，$r$为预测时间步，$n$为当前时间步，对方程离散得到

- 注意$\Phi$和$p$均采用当前已知的通量和压力。

$$
\frac{\mathbf{U}_{\mathrm{P}}^{r}-\mathbf{U}_{\mathrm{P}}^{n}}{\Delta t} V_{\mathrm{P}}+\sum \Phi_{f}^{n} \mathbf{U}_{f}^{r} - \sum \nu \nabla_{f} \mathbf{U}^{r}\left|\mathbf{S}_{f}\right| = -\sum p_{f}^{n} \mathbf{S}_{f} \tag{1}
$$

对面速度$\boldsymbol{U}_f^r$和速度的法向面梯度$\nabla_f\boldsymbol{U}^r$需要进一步离散，这里采用$Gauss linear$
$$
\mathbf{U}_{f}^{r}=\frac{\mathbf{U}_{\mathrm{P}}^{r}+\mathbf{U}_{\mathrm{N}}^{r}}{2},\quad  p_{f}^{n}=\frac{p_{\mathrm{P}}^{n}+p_{\mathrm{N}}^{n}}{2}, \quad \nabla_{f} \mathbf{U}^{r}=\frac{\mathbf{U}_{\mathrm{N}}^{r}-\mathbf{U}_{\mathrm{P}}^{r}}{|d|}
$$
至此所有未知量的面心值都用体心值表示，带入方程(1)化简后得到 **动量预测方程**
$$
A_{\mathrm{P}} \mathbf{U}_{\mathrm{P}}^{r}+\sum A_{\mathrm{N}} \mathbf{U}_{\mathrm{N}}^{r}=\boldsymbol{S}_{\mathrm{P}}^{n}-\sum \frac{p_{\mathrm{P}}^{n}+p_{\mathrm{N}}^{n}}{2} \mathbf{S}_{f}\tag{2}
$$
其中
$$
A_{\mathrm{P}}=\left(\frac{V_{p}}{\Delta t}+\sum \frac{\Phi_{f}^{n}}{2}+\sum \nu \frac{\left|\mathbf{S}_{f}\right|}{|d|}\right),\quad A_{\mathrm{N}}=\frac{\Phi_{f}^{n}}{2}-\nu \frac{\left|\mathbf{S}_{f}\right|}{|d|},\quad \boldsymbol{S}_{\mathrm{P}}^{n}=\frac{V_{p}}{\Delta t} \mathbf{U}_{\mathrm{P}}^{n} 
$$

在当前时间步中，$A_P, A_N, \boldsymbol{S}_P$都是已知量，求解 动量预测方程(2)即可获得预测速度$\mathbf{U}_{\mathrm{P}}^{r}$



### 推导压力方程

定义$H(U)$
$$
H(\boldsymbol{U}) = \boldsymbol{S}_{\mathrm{P}}-\sum A_{\mathrm{N}} \mathbf{U}_{\mathrm{N}}
$$
同样可以导出半离散方程（压力$p$不离散）：
$$
A_{P} \mathbf{U}_{P}=\mathbf{H}(\mathbf{U})-\nabla p\tag{3}
$$
即
$$
\mathbf{U}_{P}=\frac{\mathbf{H} (\mathbf{U})}{A_{P}}-\frac{1}{A_{P}} \nabla p\tag{4}
$$
由连续性方程得到**压力泊松方程**
$$
\nabla \cdot\left(\frac{\mathbf{H}(\mathbf{U})}{A_{P}}\right)  = \nabla \cdot\left(\frac{1}{A_{P}} \nabla p\right) \tag{5}
$$
对右端压力可采用Laplacian离散格式，可以求出压力$p$

为了得到通量$\Phi$， 同样在$f$面心处也有
$$
\mathbf{U}_{f}=\left(\frac{\mathbf{H} (\mathbf{U})}{A_{P}}\right)_{f}-\left(\frac{1}{A_{P}}\right)_{f}(\nabla p)_{f}
$$
所以
$$
\Phi=\mathbf{S}_f \cdot\mathbf{U}_{f}=\mathbf{S}_f \cdot\left[\left(\frac{\mathbf{H}(\mathbf{U})}{A_{P}}\right)_{f}-\left(\frac{1}{A_{P}}\right)_{f}(\nabla p)_{f}\right]\tag{6}
$$
最后得到原方程的离散形式
$$
\left\{
\begin{array}{l}
A_{P} \mathbf{U}_{P}=\mathbf{H}(\mathbf{U})-\nabla p \\
\sum \mathbf{S}_f \cdot\left[\left(\frac{1}{A_{P}}\right)_{f}(\nabla p)_{f}\right]=\sum \mathbf{S}_f \cdot\left(\frac{\mathbf{H}(\mathbf{U})}{A_{P}}\right)_{f}
\end{array}
\right. \tag{*}
$$

### PISO 算法具体流程

- #### 动量预测

  - 求解动量预测方程（2）

    ```C++
    fvVectorMatrix UEqn
        (
        fvm::ddt(U)
        + fvm::div(phi, U)
        - fvm::laplacian(nu, U)
    	);
    
    if (piso.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
    }
    ```

- #### 更新压力
  
  - 用预测的速度组建$\mathbf{H}(\mathbf{U})$
  
    ```c++
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
    surfaceScalarField phiHbyA
        (
        "phiHbyA",
        fvc::flux(HbyA)
        + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
    	);
    
    adjustPhi(phiHbyA, U, p);
    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p, U, phiHbyA, rAU);
    
    ```
  
  - 组建压力方程（5）并求解，得到更新后的压力
  
    ```c++
     while (piso.correctNonOrthogonal())
     {
         // Pressure corrector
    
         fvScalarMatrix pEqn
             (
             fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
         );
    
         pEqn.setReference(pRefCell, pRefValue);
    
         pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
    
         if (piso.finalNonOrthogonalIter())
         {
             phi = phiHbyA - pEqn.flux();
         }
     }
    ```
  
- #### 显式的速度校正
  
  - 得到新的压力后通过（6）更新通量$\Phi$
  
    ```c++
    if (piso.finalNonOrthogonalIter())
    {
    	phi = phiHbyA - pEqn.flux();
    } 
    ```
  
  - 通过（4）用新的压力显式地更新速度
  
    ```c++
    #include "continuityErrs.H"
     U = HbyA - rAU*fvc::grad(p);
     U.correctBoundaryConditions();
    ```



### 重要函数

```c++
// 解析涉及到的函数
piso.momentumPredictor()
solve
piso.correct()
volScalarField  volVectorField  fvVectorMatrix //用到哪些构造函数
fvc::flux(HbyA)
fvc::interpolate(rAU)
fvc::ddtCorr(U, phi)
adjustPhi
constrainPressure
piso.correctNonOrthogonal()
pEqn.setReference
pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())))
piso.finalNonOrthogonalIter()
U.correctBoundaryConditions()
runTime.loop() runTime.write()
```



### 一些对象和变量

```c++
runTime
piso
U p nu
UEqn UEqn.A() UEqn.H()
rAU
HbyA
phiHbyA
pEqn pEqn.flux()
```





## 附录

### icoFoam.C

```c++
#include "fvCFD.H"
#include "pisoControl.H"


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );

            adjustPhi(phiHbyA, U, p);

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAU);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );

                pEqn.setReference(pRefCell, pRefValue);

                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }
            }

            #include "continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}
```



### createFields.H

```C++
Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

```



###  Foam::fvMatrix<Type>::H() const

```C++
Foam::fvMatrix<Type>::H() const
 {
     tmp<GeometricField<Type, fvPatchField, volMesh>> tHphi
     (
         new GeometricField<Type, fvPatchField, volMesh>
         (
             IOobject
             (
                 "H("+psi_.name()+')',
                 psi_.instance(),
                 psi_.mesh(),
                 IOobject::NO_READ,
                 IOobject::NO_WRITE
             ),
             psi_.mesh(),
             dimensions_/dimVol,
             extrapolatedCalculatedFvPatchScalarField::typeName
         )
     );
     GeometricField<Type, fvPatchField, volMesh>& Hphi = tHphi.ref();
 
     // Loop over field components
     for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
     {
         scalarField psiCmpt(psi_.primitiveField().component(cmpt));
 
         scalarField boundaryDiagCmpt(psi_.size(), 0.0);
         addBoundaryDiag(boundaryDiagCmpt, cmpt);
         boundaryDiagCmpt.negate();
         addCmptAvBoundaryDiag(boundaryDiagCmpt);
 
         Hphi.primitiveFieldRef().replace(cmpt, boundaryDiagCmpt*psiCmpt);
     }
 
     Hphi.primitiveFieldRef() += lduMatrix::H(psi_.primitiveField()) + source_;
     addBoundarySource(Hphi.primitiveFieldRef());
 
     Hphi.primitiveFieldRef() /= psi_.mesh().V();
     Hphi.correctBoundaryConditions();
 
     typename Type::labelType validComponents
     (
         psi_.mesh().template validComponents<Type>()
     );
 
     for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
     {
         if (validComponents[cmpt] == -1)
         {
             Hphi.replace
             (
                 cmpt,
                 dimensionedScalar("0", Hphi.dimensions(), 0.0)
             );
         }
     }
 
     return tHphi;
 }
```



