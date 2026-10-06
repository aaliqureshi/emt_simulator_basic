function diverge39(z, m, outfile)
h = '/Users/aali27/Work/repos/powersas.m'; cd(h);
addpath(h); addpath([h '/data']); addpath([h '/util']); addpath([h '/internal']); addpath([h '/logging']);
global IS_OCTAVE Global_Settings; IS_OCTAVE = false; Global_Settings.logLevel = 'WARN';
st = runPowerSAS('pf', 'd_039_mod.m', []);
SysData = readDataFile('d_039_mod', [h '/data'], h);
for mm = [m, floor(m)]   % NR run and the HE-algebraic reference with the same integrator
  r = mytsa3(SysData, st.snapshot, [1, 0.00, 0, z, 0.5, 0.75], mm, 0.01, 1.5);
  [~, idxs] = getIndexDyn(r.SysDataBase); t = r.t; x = r.stateCurve;
  V = abs(x(idxs.vIdx, :)); d = x(idxs.deltaIdx, :); w = x(idxs.omegaIdx, :);
  fprintf('method %.1f: flag=%d msg="%s" t_end=%.3f\n', mm, r.flag, r.msg, t(end));
  for tq = [0.74 0.76 0.8 0.9 1.0 1.1 1.2 1.3 1.35 1.38 1.4 t(end)]
    k = find(t <= tq + 1e-9, 1, 'last');
    fprintf('  t=%.3f  V1=%.3f minV=%.3f  max|w|=%.4f  delta spread=%.1f deg\n', t(k), V(1,k), min(V(:,k)), max(abs(w(:,k))), (max(d(:,k)) - min(d(:,k)))*180/pi);
  end
end
end

function res = mytsa3(SysData, snapshot, faultSpec, method, dt, simLen)
simSettings = getDefaultSimSettings();
options.nlvl=15; options.taylorN=4; options.segAlpha=0.15; options.dAlpha=1; options.alphaTol=0.001;
options.diffTol=1e-6; options.diffTolMax=1e-2; options.method=method; options.simLen=simLen;
options.diffTolCtrl=1e-5; options.Efstd=1.2; options.hotStart=1;
fe = [faultSpec(:,1:end-2), 0, faultSpec(end-1); faultSpec(:,1:end-2), 1, faultSpec(end)];
simSettings.evtFault = [1 1 1; 2 2 2]; simSettings.evtFaultSpec = fe(:,1:end-1);
simSettings.eventList = [1, 0, 0, 0, 1, method, dt; 2, fe(1,end), 0, 6, 1, method, dt; ...
                         3, fe(2,end), 0, 6, 2, method, dt; 4, simLen, 0, 99, 0, method, dt];
res = runDynamicSimulationExec('d_039_mod.m', SysData, simSettings, options, snapshot);
end
