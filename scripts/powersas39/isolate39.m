function isolate39(zlist, outdir, simLen)
% Which event does Newton fail at? Methods are set per event row:
%   dyn = method of the time-domain segments, app/clr = method of the fault
%   application / clearing switch (algebraic re-init). x.0 = HE, x.1 = NR.
h = '/Users/aali27/Work/repos/powersas.m'; cd(h);
addpath(h); addpath([h '/data']); addpath([h '/util']); addpath([h '/internal']); addpath([h '/logging']);
global IS_OCTAVE Global_Settings; IS_OCTAVE = false; Global_Settings.logLevel = 'WARN';
st = runPowerSAS('pf', 'd_039_mod.m', []);
SysData = readDataFile('d_039_mod', [h '/data'], h);
combos = {'dynHE_appHE_clrNR', 1.0, 0.0, 0.1; ...
          'dynHE_appNR_clrHE', 1.0, 0.1, 0.0; ...
          'dynNR_appHE_clrHE', 1.1, 0.0, 0.0; ...
          'dynHE_appHE_clrHE', 1.0, 0.0, 0.0};
fid = fopen([outdir '/isolate39.csv'], 'w');
fprintf(fid, 'zf,combo,status,v1_pre_clear,v1_post_clear,vmin_post_clear,vmin_end\n');
for z = zlist
  for c = 1:size(combos, 1)
    name = combos{c,1}; dyn = combos{c,2}; app = floor(dyn) + combos{c,3}; clr = floor(dyn) + combos{c,4};
    status = 'error'; v1pre = NaN; v1post = NaN; vminpost = NaN; vmine = NaN;
    try
      r = mytsa2(SysData, st.snapshot, [1, 0.00, 0, z, 0.5, 0.75], dyn, app, clr, 0.01, simLen);
      [~, idxs] = getIndexDyn(r.SysDataBase); V = abs(r.stateCurve(idxs.vIdx, :)); t = r.t;
      k = find(abs(t - 0.75) < 1e-9);          % pre- and post-switch columns at t_clear
      if numel(k) >= 2; kpre = k(1); kpost = k(end); else; kpre = find(t < 0.75, 1, 'last'); kpost = kpre + 1; end
      v1pre = V(1, kpre); v1post = V(1, kpost); vminpost = min(V(:, kpost)); vmine = min(V(:, end));
      if any(~isfinite(V(:))) || t(end) < simLen - 1e-6; status = 'diverge';
      elseif vmine < 0.5; status = 'lowvoltage'; else; status = 'normal'; end
      if strcmp(name, 'dynHE_appHE_clrHE')
        xpre = r.stateCurve(:, kpre); xpost = r.stateCurve(:, kpost); tt = t;
        save(sprintf('%s/state_clear_zf%.4f.mat', outdir, z), 'xpre', 'xpost', 'idxs', 'tt', 'z');
      end
    catch ME
      status = ['error: ' strrep(ME.message, ',', ';')];
    end
    fprintf(fid, '%.4f,%s,%s,%.4f,%.4f,%.4f,%.4f\n', z, name, status, v1pre, v1post, vminpost, vmine);
    fprintf('zf=%.4f %-18s -> %s (V1 at clearing %.3f -> %.3f, vmin after %.3f, vmin end %.3f)\n', z, name, status, v1pre, v1post, vminpost, vmine);
  end
end
fclose(fid);
end

function res = mytsa2(SysData, snapshot, faultSpec, dyn, app, clr, dt, simLen)
simSettings = getDefaultSimSettings();
options.nlvl=15; options.taylorN=4; options.segAlpha=0.15; options.dAlpha=1; options.alphaTol=0.001;
options.diffTol=1e-6; options.diffTolMax=1e-2; options.method=dyn; options.simLen=simLen;
options.diffTolCtrl=1e-5; options.Efstd=1.2; options.hotStart=1;
fe = [faultSpec(:,1:end-2), 0, faultSpec(end-1); faultSpec(:,1:end-2), 1, faultSpec(end)];
simSettings.evtFault = [1 1 1; 2 2 2];
simSettings.evtFaultSpec = fe(:,1:end-1);
simSettings.eventList = [1, 0, 0, 0, 1, dyn, dt; ...
                         2, fe(1,end), 0, 6, 1, app, dt; ...
                         3, fe(2,end), 0, 6, 2, clr, dt; ...
                         4, simLen, 0, 99, 0, dyn, dt];
res = runDynamicSimulationExec('d_039_mod.m', SysData, simSettings, options, snapshot);
end
