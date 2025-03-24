PlotEmpiricalKS('maxKS','edge','parc',100:100:1000,'iFOD2')
print('EmpParcAll_A.png','-dpng')
close all

PlotEmpiricalKS('maxKS','DegCorr','parc',100:100:1000,'iFOD2')
print('EmpParcAll_B.png','-dpng')
close all
PlotEmpiricalKS('edge','DegCorr','parc',100:100:1000,'iFOD2')
print('EmpParcAll_C.png','-dpng')
close all
PlotEmpiricalKS('maxKS','edge','parc',100:100:1000,'FACT')
print('FACTEmpParcAll_A.png','-dpng')
close all
PlotEmpiricalKS('maxKS','DegCorr','parc',100:100:1000,'FACT')
print('FACTEmpParcAll_B.png','-dpng')
close all
PlotEmpiricalKS('edge','DegCorr','parc',100:100:1000,'FACT')
print('FACTEmpParcAll_C.png','-dpng')
close all