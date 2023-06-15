function out = model
%
% viga4.m
%
% Model exported on Jun 6 2023, 13:30 by COMSOL 5.6.0.280.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('D:\lucas\Estagio2023\New PINN\Studies\9 E complexo');

model.label('viga3.mph');

model.param.set('L', '0.7', 'Longueur poutre');
model.param.set('b', '2e-2', 'Largeur poutre');
model.param.set('e', '1e-2', 'Epaisseur poutre');
model.param.set('E1', '70e9', 'Module Young');
model.param.set('nu1', '0.3', 'Coefficient Poisson');
model.param.set('eta1', '0.005', 'Facteur de perte');
model.param.set('rho1', '2700', 'Masse volumique');
model.param.set('x_F', '0.1', 'Abscisse Effort ponctuel');
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
model.param.set('m', '0.32');
model.param.set('u', '0.03');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');
model.result.table.create('tbl4', 'Table');
model.result.table.create('tbl5', 'Table');
model.result.table.create('tbl6', 'Table');
model.result.table.create('tbl7', 'Table');

model.component('comp1').func.create('an1', 'Analytic');
model.component('comp1').func.create('an2', 'Analytic');
model.component('comp1').func('an1').label('E_variable');
model.component('comp1').func('an1').set('expr', '(((((1/(u*sqrt(2*pi)))*exp((-0.5)*(((x-m)/u))^2))/-40 ))+1)*E1');
model.component('comp1').func('an1').set('plotargs', {'x' '0' 'L'});
model.component('comp1').func('an2').label('Test');
model.component('comp1').func('an2').set('expr', 'E1');
model.component('comp1').func('an2').set('plotargs', {'x' '0' 'L'});

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('ls1', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls1').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls1').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls1').set('coord2', {'x_E1' '0'});
model.component('comp1').geom('geom1').create('ls2', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls2').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls2').set('coord1', {'x_E1' '0'});
model.component('comp1').geom('geom1').feature('ls2').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls2').set('coord2', {'x_E2' '0'});
model.component('comp1').geom('geom1').create('ls3', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls3').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls3').set('coord1', {'x_E2' '0'});
model.component('comp1').geom('geom1').feature('ls3').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls3').set('coord2', {'L' '0'});
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').feature('pt1').active(false);
model.component('comp1').geom('geom1').feature('pt1').set('p', {'x0' '0'});
model.component('comp1').geom('geom1').create('arr1', 'Array');
model.component('comp1').geom('geom1').feature('arr1').active(false);
model.component('comp1').geom('geom1').feature('arr1').set('fullsize', {'Nx' '1'});
model.component('comp1').geom('geom1').feature('arr1').set('displ', {'dx' '0'});
model.component('comp1').geom('geom1').feature('arr1').selection('input').set({'pt1'});
model.component('comp1').geom('geom1').create('pt2', 'Point');
model.component('comp1').geom('geom1').feature('pt2').label('Point F');
model.component('comp1').geom('geom1').feature('pt2').set('p', {'x_F' '0'});
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material.create('mat2', 'Common');
model.material.create('mat3', 'Common', '');
model.component('comp1').material('mat2').selection.set([3]);

model.component('comp1').common.create('mpf1', 'ParticipationFactors');

model.component('comp1').physics.create('beam', 'HermitianBeam', 'geom1');
model.component('comp1').physics('beam').create('fix1', 'Fixed', 0);
model.component('comp1').physics('beam').feature('fix1').selection.set([1]);
model.component('comp1').physics('beam').create('pl1', 'PointLoad', 0);
model.component('comp1').physics('beam').feature('pl1').selection.set([2]);

model.component('comp1').mesh('mesh1').create('edg2', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg2').selection.all;

model.result.table('tbl1').comments('Global Evaluation 1');
model.result.table('tbl2').comments('Global Evaluation 1');
model.result.table('tbl3').comments('Point Evaluation 1');
model.result.table('tbl4').comments('Line Maximum 1');
model.result.table('tbl5').comments('Point Evaluation 1');
model.result.table('tbl6').comments('Global Evaluation 1');
model.result.table('tbl7').comments('Point Evaluation 1');

model.component('comp1').view('view1').axis.set('xmin', -0.025015249848365784);
model.component('comp1').view('view1').axis.set('xmax', 0.5003057718276978);
model.component('comp1').view('view1').axis.set('ymin', -0.12361864745616913);
model.component('comp1').view('view1').axis.set('ymax', 0.19676025211811066);

model.component('comp1').material('mat1').propertyGroup('def').set('youngsmodulus', 'E1*(1+i*eta1)');
model.component('comp1').material('mat1').propertyGroup('def').set('poissonsratio', 'nu1');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho1');
model.component('comp1').material('mat2').active(false);
model.component('comp1').material('mat2').propertyGroup('def').set('youngsmodulus', 'E2*(1+i*eta1)');
model.component('comp1').material('mat2').propertyGroup('def').set('poissonsratio', 'nu1');
model.component('comp1').material('mat2').propertyGroup('def').set('density', 'rho1');
model.material('mat3').propertyGroup('def').set('youngsmodulus', 'an1');

model.component('comp1').physics('beam').feature('emm1').set('E', 'an1(x[1/m])');
model.component('comp1').physics('beam').feature('csd1').set('CrossSectionDefinition', 'CommonSections');
model.component('comp1').physics('beam').feature('csd1').set('hy_rect', 'e');
model.component('comp1').physics('beam').feature('csd1').set('hz_rect', 'b');
model.component('comp1').physics('beam').feature('pl1').set('Fp', {'0'; 'F'; '0'});

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
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
model.result.numerical.create('pev1', 'EvalPoint');
model.result.numerical('pev1').selection.set([5]);
model.result.numerical('pev1').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result.create('pg3', 'PlotGroup2D');
model.result.create('pg4', 'PlotGroup2D');
model.result.create('pg5', 'PlotGroup2D');
model.result.create('pg6', 'PlotGroup2D');
model.result.create('pg7', 'PlotGroup1D');
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
model.result('pg7').create('lngr1', 'LineGraph');
model.result('pg7').feature('lngr1').set('xdata', 'expr');
model.result('pg7').feature('lngr1').selection.set([1 2 3 4]);

model.nodeGroup.create('dset1beamsfgrp', 'Results');
model.nodeGroup('dset1beamsfgrp').set('type', 'plotgroup');
model.nodeGroup('dset1beamsfgrp').placeAfter('plotgroup', 'pg2');
model.nodeGroup.create('dset1beamlgrp', 'Results');
model.nodeGroup('dset1beamlgrp').set('type', 'plotgroup');
model.nodeGroup.move('dset1beamlgrp', 1);
model.nodeGroup('dset1beamlgrp').placeAfter('plotgroup', 'pg2');

model.study('std1').feature('freq').setIndex('plist', '577', 0);

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Frequency Domain');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
model.sol('sol1').feature('v1').set('cname', {'freq'});
model.sol('sol1').feature('v1').set('clist', {'577[Hz]'});
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').feature('p1').label('Parametric 1.1');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', [577]);
model.sol('sol1').feature('s1').feature('p1').set('punit', {'Hz'});
model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'auto');
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').runAll;

model.result.numerical('pev1').set('looplevelinput', {'manual'});
model.result.numerical('pev1').set('table', 'tbl7');
model.result.numerical('pev1').set('expr', {'abs(v)'});
model.result.numerical('pev1').set('unit', {'m'});
model.result.numerical('pev1').set('descr', {''});
model.result.numerical('pev1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result.numerical('pev1').setResult;
model.result('pg1').label('Mesh Plot 1');
model.result('pg1').set('inherithide', true);
model.result('pg2').label('Stress (beam)');
model.result('pg2').feature('line1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg2').feature('line1').set('linetype', 'tube');
model.result('pg2').feature('line1').set('radiusexpr', 'beam.re');
model.result('pg2').feature('line1').set('tuberadiusscaleactive', true);
model.result('pg2').feature('line1').set('colortable', 'RainbowLight');
model.result('pg2').feature('line1').set('resolution', 'normal');
model.result('pg2').feature('line1').feature('def').set('scale', 9796.953697674899);
model.result('pg2').feature('line1').feature('def').set('scaleactive', false);
model.result('pg3').label('Moment (beam)');
model.result('pg3').feature('line1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg3').feature('line1').set('linetype', 'tube');
model.result('pg3').feature('line1').set('radiusexpr', 'beam.re');
model.result('pg3').feature('line1').set('tuberadiusscaleactive', true);
model.result('pg3').feature('line1').set('colortable', 'Wave');
model.result('pg3').feature('line1').set('colortablesym', true);
model.result('pg3').feature('line1').set('resolution', 'normal');
model.result('pg4').label('Shear Force (beam)');
model.result('pg4').feature('line1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg4').feature('line1').set('linetype', 'tube');
model.result('pg4').feature('line1').set('radiusexpr', 'beam.re');
model.result('pg4').feature('line1').set('tuberadiusscaleactive', true);
model.result('pg4').feature('line1').set('colortable', 'Wave');
model.result('pg4').feature('line1').set('colortablesym', true);
model.result('pg4').feature('line1').set('resolution', 'normal');
model.result('pg5').label('Axial Force (beam)');
model.result('pg5').feature('line1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg5').feature('line1').set('linetype', 'tube');
model.result('pg5').feature('line1').set('radiusexpr', 'beam.re');
model.result('pg5').feature('line1').set('tuberadiusscaleactive', true);
model.result('pg5').feature('line1').set('colortable', 'Wave');
model.result('pg5').feature('line1').set('colortablesym', true);
model.result('pg5').feature('line1').set('resolution', 'normal');
model.result('pg6').label('Point Loads (beam)');
model.result('pg6').set('titletype', 'label');
model.result('pg6').set('frametype', 'spatial');
model.result('pg6').set('showlegendsunit', true);
model.result('pg6').feature('arpt1').label('Point Load 1');
model.result('pg6').feature('arpt1').set('expr', {'beam.pl1.F_Px' 'beam.pl1.F_Py'});
model.result('pg6').feature('arpt1').set('descr', 'Load');
model.result('pg6').feature('arpt1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg6').feature('arpt1').feature('col').set('coloring', 'gradient');
model.result('pg6').feature('arpt1').feature('col').set('topcolor', 'red');
model.result('pg6').feature('arpt1').feature('col').set('bottomcolor', 'custom');
model.result('pg6').feature('arpt1').feature('col').set('custombottomcolor', [0.5882353186607361 0.5137255191802979 0.5176470875740051]);
model.result('pg6').feature('arpt1').feature('def').set('scale', 0);
model.result('pg6').feature('arpt1').feature('def').set('scaleactive', true);
model.result('pg7').set('looplevelinput', {'manual'});
model.result('pg7').set('xlabel', 'x-coordinate (m)');
model.result('pg7').set('ylabel', 'Displacement field, y component (m)');
model.result('pg7').set('xlabelactive', false);
model.result('pg7').set('ylabelactive', false);
model.result('pg7').feature('lngr1').set('const', {'beam.refpntx' '0' 'Reference point for moment computation, x coordinate'; 'beam.refpnty' '0' 'Reference point for moment computation, y coordinate'; 'beam.refpntz' '0' 'Reference point for moment computation, z coordinate'});
model.result('pg7').feature('lngr1').set('xdataexpr', 'x');
model.result('pg7').feature('lngr1').set('xdatadescr', 'x-coordinate');
model.result('pg7').feature('lngr1').set('resolution', 'normal');

model.nodeGroup('dset1beamsfgrp').label('Section Forces (beam)');
model.nodeGroup('dset1beamsfgrp').add('plotgroup', 'pg3');
model.nodeGroup('dset1beamsfgrp').add('plotgroup', 'pg4');
model.nodeGroup('dset1beamsfgrp').add('plotgroup', 'pg5');
model.nodeGroup('dset1beamlgrp').label('Applied Loads (beam)');
model.nodeGroup('dset1beamlgrp').add('plotgroup', 'pg6');

out = model;
