function export39(zlist, outroot)
% Newton only at the clearing switch (HE elsewhere); the shadow solveAlgebraicNR
% dumps the clearing problem and Newton's result to JSON.
h = '/Users/aali27/Work/repos/powersas.m'; cd(h);
addpath(h); addpath([h '/data']); addpath([h '/util']); addpath([h '/internal']); addpath([h '/logging']);
global IS_OCTAVE Global_Settings NR_EXPORT_DIR NR_CALL; IS_OCTAVE = false; Global_Settings.logLevel = 'WARN';
st = runPowerSAS('pf', 'd_039_mod.m', []);
addpath([outroot '/shadow'], '-begin');           % shadow copy must win over internal/
disp(which('solveAlgebraicNR'));
SysData = readDataFile('d_039_mod', [h '/data'], h);
for z = zlist
  d = sprintf('%s/export_zf%.4f', outroot, z); mkdir(d);
  NR_EXPORT_DIR = d; NR_CALL = 0;
  r = mytsa2x(SysData, st.snapshot, [1, 0.00, 0, z, 0.5, 0.75], 1.0, 1.0, 1.1, 0.01, 0.8);
  [~, idxs] = getIndexDyn(r.SysDataBase); V = abs(r.stateCurve(idxs.vIdx, :)); t = r.t;
  k = find(abs(t - 0.75) < 1e-9);
  fprintf('zf=%.4f: NR calls=%d, V1 at clearing %.4f -> %.4f\n', z, NR_CALL, V(1, k(1)), V(1, k(end)));
  % HE solution of the same clearing switch, for comparison
  NR_EXPORT_DIR = [];
  r2 = mytsa2x(SysData, st.snapshot, [1, 0.00, 0, z, 0.5, 0.75], 1.0, 1.0, 1.0, 0.01, 0.8);
  [~, idxs2] = getIndexDyn(r2.SysDataBase); t2 = r2.t; k2 = find(abs(t2 - 0.75) < 1e-9);
  Vhe = r2.stateCurve(idxs2.vIdx, k2(end)); Vpre = r2.stateCurve(idxs2.vIdx, k2(1));
  fid = fopen([d '/he_clear.json'], 'w');
  fprintf(fid, '%s', jsonencode(struct('Vr', real(Vhe), 'Vi', imag(Vhe), 'Vpre_r', real(Vpre), 'Vpre_i', imag(Vpre)))); fclose(fid);
end
end

function res = mytsa2x(SysData, snapshot, faultSpec, dyn, app, clr, dt, simLen)
simSettings = getDefaultSimSettings();
options.nlvl=15; options.taylorN=4; options.segAlpha=0.15; options.dAlpha=1; options.alphaTol=0.001;
options.diffTol=1e-6; options.diffTolMax=1e-2; options.method=dyn; options.simLen=simLen;
options.diffTolCtrl=1e-5; options.Efstd=1.2; options.hotStart=1;
fe = [faultSpec(:,1:end-2), 0, faultSpec(end-1); faultSpec(:,1:end-2), 1, faultSpec(end)];
simSettings.evtFault = [1 1 1; 2 2 2];
simSettings.evtFaultSpec = fe(:,1:end-1);
simSettings.eventList = [1, 0, 0, 0, 1, dyn, dt; 2, fe(1,end), 0, 6, 1, app, dt; ...
                         3, fe(2,end), 0, 6, 2, clr, dt; 4, simLen, 0, 99, 0, dyn, dt];
res = runDynamicSimulationExec('d_039_mod.m', SysData, simSettings, options, snapshot);
end
