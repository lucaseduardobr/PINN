function out = model
%
% Viga.m
%
% Model exported on Apr 17 2023, 14:33 by COMSOL 5.6.0.280.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('G:\Utilisateurs\lcoelho2\Desktop\Estagio rafael');

model.label('Poutre_Mael.mph');

data.L = 0.7;
data.rho = 2700;
data.b = 2e-2;
data.e = 1e-2;
data.e2 = 0.8*1e-2;
data.S = data.b*data.e;
data.nu = 0.3;
data.eta1 = 0;
data.I = data.b*data.e^3/(12); %*(1-nu)
data.E1 = 70e9;
data.E2 = 0.8*70e9; %0.8*  ;1*data.E1
data.F = 1;
data.x_F=0.1;
data.freq_min = 500;
data.freq_max = 600;




model.param.set('L', '0.7', 'Longueur poutre');
model.param.set('b', '2e-2', 'Largeur poutre');
model.param.set('e', '1e-2', 'Epaisseur poutre');
model.param.set('E1', '70e9', 'Module Young');
model.param.set('nu1', '0.3', 'Coefficient Poisson');
model.param.set('eta1', '0', 'Facteur de perte');
model.param.set('rho1', '2700', 'Masse volumique');
model.param.set('x_F', '0.6', 'Abscisse Effort ponctuel');
model.param.set('F', '1', 'Amplitude effort');
model.param.set('x_E1', '0.4', 'Abscisse changement module');
model.param.set('x_E2', '0.45', 'Abscisse fin changement module');
model.param.set('E2', '0.8*E1', ['Module Young modifi' native2unicode(hex2dec({'00' 'e9'}), 'unicode') ]);
model.param.set('x0', '0.3', 'Abscisse premier point mesure');
model.param.set('dx', '0.031', 'Pas de mesure');
model.param.set('Nx', '10', 'Nb points de mesure');
model.param.set('freq_min', '500', ['Fr' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quence min']);
model.param.set('freq_max', '600', ['Fr' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quence max']);
model.param.set('freq_step', '2', ['Pas en fr' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quence']);
model.param.label('Parameters 1');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').label('Component 1');

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').label('Geometry 1');
model.component('comp1').geom('geom1').create('ls1', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls1').label('Line Segment 1');
model.component('comp1').geom('geom1').feature('ls1').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls1').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls1').set('coord2', {'x_E1' '0'});
model.component('comp1').geom('geom1').create('ls2', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls2').label('Line Segment 2');
model.component('comp1').geom('geom1').feature('ls2').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls2').set('coord1', {'x_E1' '0'});
model.component('comp1').geom('geom1').feature('ls2').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls2').set('coord2', {'x_E2' '0'});
model.component('comp1').geom('geom1').create('ls3', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls3').label('Line Segment 3');
model.component('comp1').geom('geom1').feature('ls3').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls3').set('coord1', {'x_E2' '0'});
model.component('comp1').geom('geom1').feature('ls3').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls3').set('coord2', {'L' '0'});
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').feature('pt1').active(false);
model.component('comp1').geom('geom1').feature('pt1').set('p', {'x0' '0'});
model.component('comp1').geom('geom1').create('arr1', 'Array');
model.component('comp1').geom('geom1').feature('arr1').active(false);
model.component('comp1').geom('geom1').feature('arr1').label('Array 1');
model.component('comp1').geom('geom1').feature('arr1').set('fullsize', {'Nx' '1'});
model.component('comp1').geom('geom1').feature('arr1').set('displ', {'dx' '0'});
model.component('comp1').geom('geom1').feature('arr1').selection('input').set({'pt1'});
model.component('comp1').geom('geom1').create('pt2', 'Point');
model.component('comp1').geom('geom1').feature('pt2').label('Point F');
model.component('comp1').geom('geom1').feature('pt2').set('p', {'x_F' '0'});
model.component('comp1').geom('geom1').feature('fin').label('Form Union');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat1').selection.set([1 3 4]);
model.component('comp1').material('mat2').selection.set([2]);

model.component('comp1').physics.create('beam', 'HermitianBeam', 'geom1');
model.component('comp1').physics('beam').create('fix1', 'Fixed', 0);
model.component('comp1').physics('beam').feature('fix1').selection.set([1]);
model.component('comp1').physics('beam').create('pl1', 'PointLoad', 0);
model.component('comp1').physics('beam').feature('pl1').selection.set([4]);

model.component('comp1').mesh('mesh1').create('edg2', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg2').selection.all;

model.thermodynamics.label('Thermodynamics');

model.component('comp1').view('view1').label('View 1');
model.component('comp1').view('view1').axis.label('Axis');
model.component('comp1').view('view1').axis.set('xmin', -0.23215952515602112);
model.component('comp1').view('view1').axis.set('xmax', 0.9321595430374146);
model.component('comp1').view('view1').axis.set('ymin', -0.3937753140926361);
model.component('comp1').view('view1').axis.set('ymax', 0.3937753140926361);

model.component('comp1').material('mat1').label('Material 1');
model.component('comp1').material('mat1').propertyGroup('def').label('Basic');
model.component('comp1').material('mat1').propertyGroup('def').set('youngsmodulus', 'E1*(1+i*eta1)');
model.component('comp1').material('mat1').propertyGroup('def').set('poissonsratio', 'nu1');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho1');
model.component('comp1').material('mat2').label('Material 2');
model.component('comp1').material('mat2').propertyGroup('def').label('Basic');
model.component('comp1').material('mat2').propertyGroup('def').set('youngsmodulus', 'E2*(1+i*eta1)');
model.component('comp1').material('mat2').propertyGroup('def').set('poissonsratio', 'nu1');
model.component('comp1').material('mat2').propertyGroup('def').set('density', 'rho1');

model.component('comp1').coordSystem('sys1').label('Boundary System 1');

model.common('cminpt').label('Default Model Inputs');

model.component('comp1').physics('beam').label('Beam');
model.component('comp1').physics('beam').feature('emm1').label('Linear Elastic Material 1');
model.component('comp1').physics('beam').feature('emm1').featureInfo('info').label('Equation View');
model.component('comp1').physics('beam').feature('csd1').set('CrossSectionDefinition', 'CommonSections');
model.component('comp1').physics('beam').feature('csd1').set('hy_rect', 'e');
model.component('comp1').physics('beam').feature('csd1').set('hz_rect', 'b');
model.component('comp1').physics('beam').feature('csd1').label('Cross-Section Data 1');
model.component('comp1').physics('beam').feature('csd1').featureInfo('info').label('Equation View');
model.component('comp1').physics('beam').feature('free1').label('Free 1');
model.component('comp1').physics('beam').feature('free1').featureInfo('info').label('Equation View');
model.component('comp1').physics('beam').feature('init1').label('Initial Values 1');
model.component('comp1').physics('beam').feature('init1').featureInfo('info').label('Equation View');
model.component('comp1').physics('beam').feature('fix1').label('Fixed Constraint 1');
model.component('comp1').physics('beam').feature('fix1').featureInfo('info').label('Equation View');
model.component('comp1').physics('beam').feature('pl1').set('Fp', {'0'; 'F'; '0'});
model.component('comp1').physics('beam').feature('pl1').label('Point Load 1');
model.component('comp1').physics('beam').feature('pl1').featureInfo('info').label('Equation View');

model.component('comp1').mesh('mesh1').label('Mesh 1');
model.component('comp1').mesh('mesh1').feature('size').label('Size');
model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('edg2').label('Edge 2');
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('freq', 'Frequency');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.dataset.create('mesh1', 'Mesh');
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result.create('pg3', 'PlotGroup2D');
model.result.create('pg4', 'PlotGroup2D');
model.result.create('pg5', 'PlotGroup2D');
model.result.create('pg6', 'PlotGroup2D');
model.result('pg1').set('data', 'mesh1');
model.result('pg2').create('line1', 'Line');
model.result('pg2').feature('line1').set('expr', 'beam.mises');
model.result('pg2').feature('line1').create('def', 'Deform');
model.result('pg3').create('line1', 'Line');
model.result('pg3').feature('line1').set('expr', 'beam.Mzl');
model.result('pg3').feature('line1').create('def', 'Deform');
model.result('pg4').create('line1', 'Line');
model.result('pg4').feature('line1').set('expr', 'beam.Tyl');
model.result('pg4').feature('line1').create('def', 'Deform');
model.result('pg5').create('line1', 'Line');
model.result('pg5').feature('line1').set('expr', 'beam.Nxl');
model.result('pg5').feature('line1').create('def', 'Deform');
model.result('pg6').create('arpt1', 'ArrowPoint');
model.result('pg6').feature('arpt1').create('col', 'Color');
model.result('pg6').feature('arpt1').create('def', 'Deform');
model.result('pg6').feature('arpt1').feature('col').set('expr', 'beam.pl1.F_P_Mag');

model.nodeGroup.create('dset1beamsfgrp', 'Results');
model.nodeGroup('dset1beamsfgrp').set('type', 'plotgroup');
model.nodeGroup('dset1beamsfgrp').placeAfter('plotgroup', 'pg2');
model.nodeGroup.create('dset1beamlgrp', 'Results');
model.nodeGroup('dset1beamlgrp').set('type', 'plotgroup');
model.nodeGroup.move('dset1beamlgrp', 1);
model.nodeGroup('dset1beamlgrp').placeAfter('plotgroup', 'pg2');

model.study('std1').label('Study 1');
model.study('std1').feature('freq').label('Frequency Domain');
model.study('std1').feature('freq').set('plist', 'range(freq_min,freq_step,freq_max)');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label(['Compilation des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations: Frequency Domain']);
model.sol('sol1').feature('v1').label('Dependent Variables 1');
model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
model.sol('sol1').feature('v1').set('cname', {'freq'});
model.sol('sol1').feature('v1').set('clist', {'range(freq_min,freq_step,freq_max)[Hz]'});
model.sol('sol1').feature('s1').label('Stationary Solver 1');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').feature('p1').label('Parametric 1');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'range(freq_min,freq_step,freq_max)'});
model.sol('sol1').feature('s1').feature('p1').set('punit', {'Hz'});
model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'auto');
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1');
model.sol('sol1').runAll;

model.result.label('Results');
model.result.dataset('mesh1').label('Mesh 1');
model.result('pg1').label('Mesh Plot 1');
model.result('pg1').set('inherithide', true);
model.result('pg2').label('Stress (beam)');
model.result('pg2').feature('line1').label('Line 1');
model.result('pg2').feature('line1').set('descr', 'von Mises stress');
model.result('pg2').feature('line1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg2').feature('line1').set('linetype', 'tube');
model.result('pg2').feature('line1').set('radiusexpr', 'beam.re');
model.result('pg2').feature('line1').set('tuberadiusscaleactive', true);
model.result('pg2').feature('line1').set('colortable', 'RainbowLight');
model.result('pg2').feature('line1').set('resolution', 'normal');
model.result('pg2').feature('line1').feature('def').label('Deformation');
model.result('pg2').feature('line1').feature('def').set('descr', 'Displacement field');
model.result('pg2').feature('line1').feature('def').set('scale', 12236.963710231379);
model.result('pg2').feature('line1').feature('def').set('scaleactive', false);
model.result('pg3').label('Moment (beam)');
model.result('pg3').feature('line1').label('Line 1');
model.result('pg3').feature('line1').set('descr', 'Bending moment, local z direction');
model.result('pg3').feature('line1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg3').feature('line1').set('linetype', 'tube');
model.result('pg3').feature('line1').set('radiusexpr', 'beam.re');
model.result('pg3').feature('line1').set('tuberadiusscaleactive', true);
model.result('pg3').feature('line1').set('colortable', 'Wave');
model.result('pg3').feature('line1').set('colortablesym', true);
model.result('pg3').feature('line1').set('resolution', 'normal');
model.result('pg3').feature('line1').feature('def').label('Deformation');
model.result('pg3').feature('line1').feature('def').set('descr', 'Displacement field');
model.result('pg4').label('Shear Force (beam)');
model.result('pg4').feature('line1').label('Line 1');
model.result('pg4').feature('line1').set('descr', 'Shear force, local y direction');
model.result('pg4').feature('line1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg4').feature('line1').set('linetype', 'tube');
model.result('pg4').feature('line1').set('radiusexpr', 'beam.re');
model.result('pg4').feature('line1').set('tuberadiusscaleactive', true);
model.result('pg4').feature('line1').set('colortable', 'Wave');
model.result('pg4').feature('line1').set('colortablesym', true);
model.result('pg4').feature('line1').set('resolution', 'normal');
model.result('pg4').feature('line1').feature('def').label('Deformation');
model.result('pg4').feature('line1').feature('def').set('descr', 'Displacement field');
model.result('pg5').label('Axial Force (beam)');
model.result('pg5').feature('line1').label('Line 1');
model.result('pg5').feature('line1').set('descr', 'Local axial force');
model.result('pg5').feature('line1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg5').feature('line1').set('linetype', 'tube');
model.result('pg5').feature('line1').set('radiusexpr', 'beam.re');
model.result('pg5').feature('line1').set('tuberadiusscaleactive', true);
model.result('pg5').feature('line1').set('colortable', 'Wave');
model.result('pg5').feature('line1').set('colortablesym', true);
model.result('pg5').feature('line1').set('resolution', 'normal');
model.result('pg5').feature('line1').feature('def').label('Deformation');
model.result('pg5').feature('line1').feature('def').set('descr', 'Displacement field');
model.result('pg6').label('Point Loads (beam)');
model.result('pg6').set('titletype', 'label');
model.result('pg6').set('frametype', 'spatial');
model.result('pg6').set('showlegendsunit', true);
model.result('pg6').feature('arpt1').label('Point Load 1');
model.result('pg6').feature('arpt1').set('expr', {'beam.pl1.F_Px' 'beam.pl1.F_Py'});
model.result('pg6').feature('arpt1').set('descr', 'Load');
model.result('pg6').feature('arpt1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg6').feature('arpt1').feature('col').label('Color Expression');
model.result('pg6').feature('arpt1').feature('col').set('descr', 'Load magnitude');
model.result('pg6').feature('arpt1').feature('col').set('coloring', 'gradient');
model.result('pg6').feature('arpt1').feature('col').set('topcolor', 'red');
model.result('pg6').feature('arpt1').feature('col').set('bottomcolor', 'custom');
model.result('pg6').feature('arpt1').feature('col').set('custombottomcolor', [0.5882353186607361 0.5137255191802979 0.5176470875740051]);
model.result('pg6').feature('arpt1').feature('def').label('Deformation');
model.result('pg6').feature('arpt1').feature('def').set('descr', 'Displacement field');
model.result('pg6').feature('arpt1').feature('def').set('scale', 0);
model.result('pg6').feature('arpt1').feature('def').set('scaleactive', true);

model.nodeGroup('dset1beamsfgrp').label('Section Forces (beam)');
model.nodeGroup('dset1beamsfgrp').add('plotgroup', 'pg3');
model.nodeGroup('dset1beamsfgrp').add('plotgroup', 'pg4');
model.nodeGroup('dset1beamsfgrp').add('plotgroup', 'pg5');
model.nodeGroup('dset1beamlgrp').label('Applied Loads (beam)');
model.nodeGroup('dset1beamlgrp').add('plotgroup', 'pg6');

model.label('Poutre_Mael.mph');

model.component('comp1').geom('geom1').run('pt2');
model.component('comp1').geom('geom1').create('pt3', 'Point');
model.component('comp1').geom('geom1').feature('pt3').set('p', [0.09 0]);
model.component('comp1').geom('geom1').run('pt3');
model.component('comp1').geom('geom1').feature('pt3').setIndex('p', 0.1, 0);
model.component('comp1').geom('geom1').runPre('fin');
model.component('comp1').geom('geom1').run('fin');
model.component('comp1').geom('geom1').feature.remove('pt3');

model.param.set('x_F', '0.1');

model.component('comp1').geom('geom1').run('fin');

model.component('comp1').material('mat2').active(false);
model.component('comp1').material('mat1').selection.set([1 2 3 4]);

model.sol('sol1').study('std1');

model.study('std1').feature('freq').set('notlistsolnum', 1);
model.study('std1').feature('freq').set('notsolnum', '1');
model.study('std1').feature('freq').set('listsolnum', 1);
model.study('std1').feature('freq').set('solnum', '1');

model.sol('sol1').feature.remove('s1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'freq');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'freq');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').feature.remove('pDef');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'range(freq_min,freq_step,freq_max)'});
model.sol('sol1').feature('s1').feature('p1').set('punit', {'Hz'});
model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'auto');
model.sol('sol1').feature('s1').feature('p1').set('pdistrib', 'off');
model.sol('sol1').feature('s1').feature('p1').set('plot', 'off');
model.sol('sol1').feature('s1').feature('p1').set('plotgroup', 'pg1');
model.sol('sol1').feature('s1').feature('p1').set('probesel', 'all');
model.sol('sol1').feature('s1').feature('p1').set('probes', {});
model.sol('sol1').feature('s1').feature('p1').set('control', 'freq');
model.sol('sol1').feature('s1').set('control', 'freq');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg2').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 1, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 2, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 3, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 4, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 5, 0);
model.result('pg2').run;
model.result('pg1').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 1, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 2, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 3, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 4, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 5, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 6, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 7, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 8, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 9, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 10, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 11, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 12, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 13, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 14, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 15, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 16, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 17, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 18, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 19, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 20, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 21, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 22, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 23, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 24, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 25, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 26, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 27, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 28, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 29, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 30, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 31, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 32, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 33, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 34, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 35, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 51, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 50, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 49, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 48, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 47, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 46, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 45, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 44, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 43, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 42, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 41, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 40, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 39, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 38, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 37, 0);
model.result('pg2').run;
model.result('pg2').setIndex('looplevel', 36, 0);
model.result('pg2').run;
model.result('pg2').run;

%getting only a part of the beam
%step dx start at 'x0'
coord=str2double(model.param.get('x0')):str2double(model.param.get('dx')):(str2double(model.param.get('x0'))+(str2double(model.param.get('Nx'))-1)*str2double(model.param.get('dx')));
coord=[coord;zeros(size(coord))];
%edim is the dimention 1D 2D 3D
disp = mphinterp(model,'v','edim',1,'coord',coord,'solnum',1); % ,'solnum',1 pour sortir uniquement la 1e fréquence or the number of the solution
%getting the whole beam
coord_tot = 0:1e-2*str2double(model.param.get('dx')):str2double(model.param.get('L'));
coord_tot=[coord_tot;zeros(size(coord_tot))];
%coord is the points you wanna know the value
disp_tot = mphinterp(model,'v','edim',1,'coord',coord_tot,'solnum',1);

figure
plot(coord(1,:),disp(1,:))
X_COMSOL = coord(1,:);
W_COMSOL = disp(1,:);
save("displacement_part","X_COMSOL","W_COMSOL","data");

figure
plot(coord_tot(1,:),disp_tot(1,:))
X_COMSOL_tot = coord_tot(1,:);
W_COMSOL_tot = disp_tot(1,:);
save("displacement_tot","X_COMSOL_tot","W_COMSOL_tot","data");

model.label('Viga.mph');

out = model;
