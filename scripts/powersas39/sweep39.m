function sweep39(zlist, mlist, outfile, simLen)
h = '/Users/aali27/Work/repos/powersas.m'; cd(h);
addpath(h); addpath([h '/data']); addpath([h '/util']); addpath([h '/internal']); addpath([h '/logging']);
global IS_OCTAVE Global_Settings; IS_OCTAVE = false; Global_Settings.logLevel = 'WARN';
st = runPowerSAS('pf', 'd_039_mod.m', []);
SysData = readDataFile('d_039_mod', [h '/data'], h);
fid = fopen(outfile, 'w');
fprintf(fid, 'zf,method,status,t_end,v1_fault,v1_after,v1_end,vmin_end,n_low_end,msg\n');
for z = zlist
  for m = mlist
    fs = [1, 0.00, 0, z, 0.5, 0.75];
    msg = ''; t_end = NaN; v1f = NaN; v1a = NaN; v1e = NaN; vmin = NaN; nlow = NaN; status = 'error';
    try
      r = mytsa(SysData, st.snapshot, fs, m, 0.01, simLen);
      [~, idxs] = getIndexDyn(r.SysDataBase);
      V = abs(r.stateCurve(idxs.vIdx, :)); t = r.t;
      t_end = t(end);
      v1f = V(1, find(t <= 0.745, 1, 'last')); v1a = V(1, find(t >= 0.76, 1, 'first'));
      if isempty(v1a); v1a = NaN; end
      v1e = V(1, end); vmin = min(V(:, end)); nlow = sum(V(:, end) < 0.5);
      if any(~isfinite(V(:))) || t_end < simLen - 1e-6
        status = 'diverge';
      elseif vmin < 0.5
        status = 'lowvoltage';
      else
        status = 'normal';
      end
    catch ME
      msg = strrep(ME.message, ',', ';'); msg = strrep(msg, newline, ' ');
    end
    fprintf(fid, '%.4f,%.1f,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%s\n', z, m, status, t_end, v1f, v1a, v1e, vmin, nlow, msg);
    fprintf('zf=%.4f m=%.1f -> %s (t_end=%.3f v1: %.3f -> %.3f, end %.3f, vmin %.3f)\n', z, m, status, t_end, v1f, v1a, v1e, vmin);
  end
end
fclose(fid);
end

function res = mytsa(SysData, snapshot, faultSpec, method, dt, simLen)
% Copy of calTSA with an explicit step size for the numerical-integration methods.
simSettings = getDefaultSimSettings();
options.nlvl=15; options.taylorN=4; options.segAlpha=0.15; options.dAlpha=1; options.alphaTol=0.001;
options.diffTol=1e-6; options.diffTolMax=1e-2; options.method=method; options.simLen=simLen;
options.diffTolCtrl=1e-5; options.Efstd=1.2; options.hotStart=1;
faultAdd = [faultSpec(:,1:end-2), zeros(size(faultSpec,1),1), faultSpec(:,end-1)];
faultClear = [faultSpec(:,1:end-2), ones(size(faultSpec,1),1), faultSpec(:,end)];
fe = [faultAdd; faultClear]; [~, is] = sort(fe(:,end)); fe = fe(is,:);
ft = [0; fe(:,end)]; dft = ft(2:end) - ft(1:end-1); agg = [find(dft ~= 0); size(ft,1)]; n = size(agg,1) - 1;
simSettings.evtFault = [(1:n)', agg(1:n), agg(2:n+1)-1];
simSettings.evtFaultSpec = fe(:,1:end-1);
simSettings.eventList(1,6) = method; simSettings.eventList(1,7) = dt;
simSettings.eventList = [simSettings.eventList; [(2:n+1)', fe(agg(1:n),end), zeros(n,1), 6*ones(n,1), (1:n)', repmat(simSettings.eventList(1,[6,7]), n, 1)]];
simSettings.eventList = [simSettings.eventList; [n+2, simLen, 0, 99, 0, simSettings.eventList(1,[6,7])]];
res = runDynamicSimulationExec('d_039_mod.m', SysData, simSettings, options, snapshot);
end
