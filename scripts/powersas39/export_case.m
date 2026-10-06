function export_case(dataName, lineIdx, faultBusOrig, z, tOn, tOff, simLen, outdir)
% Export the fault-application and fault-clearing re-init problems of a PowerSAS
% case. Three runs: HE everywhere (reference states), NR only at application,
% NR only at clearing. The shadow solveAlgebraicNR dumps each NR problem to JSON.
h = '/Users/aali27/Work/repos/powersas.m'; cd(h);
addpath(h); addpath([h '/data']); addpath([h '/util']); addpath([h '/internal']); addpath([h '/logging']);
global IS_OCTAVE Global_Settings NR_EXPORT_DIR NR_CALL; IS_OCTAVE = false; Global_Settings.logLevel = 'WARN';
tic; st = runPowerSAS('pf', [dataName '.m'], []); fprintf('PF: %.1f s\n', toc);
addpath(fileparts(mfilename('fullpath')), '-begin'); addpath([fileparts(mfilename('fullpath')) '/shadow'], '-begin');
fprintf('solver: %s\n', which('solveAlgebraicNR'));
SysData = readDataFile(dataName, [h '/data'], h);
fb = find(SysData.bus(:,1) == faultBusOrig);
fprintf('line %d: %d-%d, fault bus index %d\n', lineIdx, SysData.line(lineIdx,1), SysData.line(lineIdx,2), fb);
mkdir(outdir);
fs = [lineIdx, 0.00, 0, z, tOn, tOff];
runs = {'he', 0.0, 0.0, 0.0; 'appNR', 0.0, 0.1, 0.0; 'clrNR', 0.0, 0.0, 0.1};
for c = 1:size(runs, 1)
  name = runs{c,1}; d = [outdir '/' name]; mkdir(d);
  NR_EXPORT_DIR = d; NR_CALL = 0; tic;
  try
    r = mytsa4(dataName, SysData, st.snapshot, fs, runs{c,2}, runs{c,3}, runs{c,4}, 0.01, simLen);
    [~, idxs] = getIndexDyn(r.SysDataBase); t = r.t; V = r.stateCurve(idxs.vIdx, :);
    ka = find(abs(t - tOn) < 1e-9); kc = find(abs(t - tOff) < 1e-9);
    out.t = t; out.flag = r.flag; out.msg = r.msg; out.fb = fb;
    out.Va_pre_r = real(V(:, ka(1))); out.Va_pre_i = imag(V(:, ka(1))); out.Va_post_r = real(V(:, ka(end))); out.Va_post_i = imag(V(:, ka(end)));
    if ~isempty(kc)
      out.Vc_pre_r = real(V(:, kc(1))); out.Vc_pre_i = imag(V(:, kc(1))); out.Vc_post_r = real(V(:, kc(end))); out.Vc_post_i = imag(V(:, kc(end)));
    end
    fid = fopen([d '/states.json'], 'w'); fprintf(fid, '%s', jsonencode(out)); fclose(fid);
    Vm = abs(V);
    fprintf('%-6s flag=%d t_end=%.3f (%.0f s) | Vfb app %.3f->%.3f', name, r.flag, t(end), toc, Vm(fb, ka(1)), Vm(fb, ka(end)));
    if ~isempty(kc); fprintf(' | clr %.3f->%.3f minV %.3f', Vm(fb, kc(1)), Vm(fb, kc(end)), min(Vm(:, kc(end)))); end
    fprintf(' | msg: %s\n', r.msg);
  catch ME
    fprintf('%-6s ERROR (%.0f s): %s\n', name, toc, ME.message);
  end
end
end

function res = mytsa4(dataName, SysData, snapshot, faultSpec, dyn, app, clr, dt, simLen)
simSettings = getDefaultSimSettings();
options.nlvl=15; options.taylorN=4; options.segAlpha=0.15; options.dAlpha=1; options.alphaTol=0.001;
options.diffTol=1e-6; options.diffTolMax=1e-2; options.method=dyn; options.simLen=simLen;
options.diffTolCtrl=1e-5; options.Efstd=1.2; options.hotStart=1;
fe = [faultSpec(:,1:end-2), 0, faultSpec(end-1); faultSpec(:,1:end-2), 1, faultSpec(end)];
simSettings.evtFault = [1 1 1; 2 2 2]; simSettings.evtFaultSpec = fe(:,1:end-1);
simSettings.eventList = [1, 0, 0, 0, 1, dyn, dt; 2, fe(1,end), 0, 6, 1, app, dt; ...
                         3, fe(2,end), 0, 6, 2, clr, dt; 4, simLen, 0, 99, 0, dyn, dt];
res = runDynamicSimulationExec([dataName '.m'], SysData, simSettings, options, snapshot);
end
