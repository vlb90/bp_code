function [uZverevBP, y, sci] = noiseMetOp(data, mode, gg)

    close all;
    
    % Directories
    
    figFolder = fullfile(data.dirSave, 'figs'); %mac
    
    %figFolder = 'C:\Users\vlb\Dropbox\Doutorado\Research\backpropagation\carrano\bp_w_noise\figs\noise\'; %win

    srcPath = strsplit(data.simFile, '\');
    
    figName = strrep(srcPath{end}, 'sim_', 'noise_');
    
    % Model SNR (including antenna pattern attenuation)
    % Riccardo's presentation (mid-latitude)
    
    slta_ref = [10 100 200 300 400 500 600].*1e3;
 
    snr_ref = [710 700 670 590 500 400 290].*1; % original in manuscript

%     snr_ref = 1000.*ones(1,length(slta_ref)); % flat SNR (without antenna radiation issue)
    
    fsWOP = round(3.2e3/mean(diff(data.Y))); % Sampling frequency on PS
    
    fsMeas  = 1;                             % Sampling freq. receiver/data   
    
    % GRAS_1B_M02_20200801020319Z_20200801020926Z_N_T_20210507122630Z_G29_NN.nc
    
    if strcmp(mode, 'meas')
    
        %load('C:\Users\vlb\Dropbox\Doutorado\Research\backpropagation\carrano\ref_GRAS_1B_M02_20200801020319Z_20200801020926Z_N_T_20210507122630Z_G29_NN.mat'); %win

        load('/Users/vlb/Dropbox/Doutorado/Research/backpropagation/carrano/ref_GRAS_1B_M02_20200801020319Z_20200801020926Z_N_T_20210507122630Z_G29_NN.mat');
    
        slta_ref = slta;
    
        snr_ref = snr_1c;
        
        clear slta snr_1c;
        
    else
 
        [fSNR, ~] = polyfit(slta_ref./1e3, snr_ref, 4);

        figure(170); plot(polyval(fSNR, slta_ref./1e3), slta_ref./1e3, 'k', 'LineWidth', 2); hold on
        
        plot(snr_ref, slta_ref./1e3, 'b--' , 'LineWidth', 1.2);
    
        xlim([-20 820]);
    
        xticks(0:100:800)
    
        ylim([35 600]);
    
        yticks([35 100:100:600]);
    
        ylabel('SLTA (km)');
    
        xlabel('SNR (V/V)');
    
        grid on;
    
        set(gcf, 'Position', [445.4000  426.2000  314.8000  420.0000]);
    
        set(gca, 'FontName', 'Courier');
    
        close(170);
        
    end
    
    % ---------------------------------------------------------------------
    
    %SLTA @ last screen
    
    z = data.Z(end-3);
    
    m = ((data.Y.*data.satellitePositionZ-data.satellitePositionY.*z)./(data.satellitePositionZ-z));
    
    k = (data.Y-data.satellitePositionY)./(z-data.satellitePositionZ);

    slta_full = m./(sqrt(1+k.^2))-data.Re;
    
    ix = 1:length(slta_full); %slta_full >= 80e3 & slta_full<=920e3;
    
    slta = slta_full(ix);
    
    y = data.Y(ix);
    
    u_signal = data.uScreenSampled(ix);
    
    % ---------------------------------------------------------------------
    
    % Filtering
    
    nPts = round(10e3/mean(diff(y)));  % Filter length: 10 km (outer scale)
    
    nPts_snr = round(20e3/mean(diff(slta_ref)));
    
    filtered_u_sq = u_signal.*conj(u_signal);
    
    snr_filtered = snr_ref;
    
    figure; plot(slta, filtered_u_sq, 'b', 'LineWidth', 1.2); hold on
    
    if strcmp(mode, 'meas') % fix added on 22-08-08 to run the code with generic SNR
    
        for kk = 1:3 % 3 pass, 2nd order Savitsky-Golay filter

            snr_filtered = sgolayfilt(snr_filtered, 2, nPts_snr);    % Reproduce figure in Appendix AMT

            filtered_u_sq = sgolayfilt(filtered_u_sq(:), 2, nPts);
    
        end
        
    else
        
        snr_filtered = snr_ref;
        
        for kk = 1:3
        
            filtered_u_sq = sgolayfilt(filtered_u_sq(:), 2, nPts);
            
        end
            
    end
    
    plot(slta, filtered_u_sq, 'r--');
    
    close all;
    
    [~, ih] = max(snr_filtered);
    
    iBelow = slta./1e3 <= slta_ref(ih)./1e3;
    
    iEnd = find(slta < max(slta_ref), 1, 'last');
    
    iAbove = slta(1:iEnd)./1e3 > slta_ref(ih)./1e3;
    
    iExtra = slta./1e3 > max(slta_ref)./1e3;
    
    [~, ihMax] = max(slta_ref);

    snr_interpolated(iBelow) = max(snr_filtered); % snr_ref

    snr_interpolated(iAbove) = interp1(slta_ref, snr_filtered, slta(iAbove)); % snr_ref
    
    snr_interpolated(iExtra) = snr_filtered(ihMax); % snr_ref
    
    sigma_0 = sqrt((abs(filtered_u_sq(:))./(snr_interpolated(:).^2))*(fsWOP/fsMeas));
    
% % % %     sigma_0 = zeros(length(snr_interpolated),1);

    if ~isnan(gg)
        
        rng(gg);
        
    end
    
    u_noise = (sigma_0./sqrt(2)).*(randn(length(y),1)+1i*randn(length(y),1));
    
    i_slta = slta./1e3 >= 280 & slta./1e3 <= 355;
    
    pd = fitdist(abs(sigma_0(i_slta)), 'Normal');

    [ci] = paramci(pd, 'Alpha', 0.05);
    
    mu_noise = pd.mu;
    
    noiseLower = ci(1,1);
    
    noiseUpper = ci(2,1);
    
    slta_int = 40:5:600;
    
    sigma_0_int = spline(slta./1e3, abs(sigma_0), slta_int);
    
    uZverevBP = (u_signal(:) + u_noise);
    
    % Noise amplitude figure
    
    figure; plot(sigma_0, slta./1e3, 'k');
    ylim([40 600])
    yticks(100:100:600)
    xlabel('\sigma_0');
    ylabel('SLTA (km)');
    set(gcf, 'Position', [445.4000  426.2000  314.8000  420.0000]);
    set(gca, 'FontName', 'Courier');
    xticks(0.04:0.01:0.08)
    set(gca, 'XGrid', 'on')
    set(gca, 'XMinorTick', 'on')
    set(gca, 'XMinorGrid', 'on')
    set(gca, 'YMinorGrid', 'on')

    % ---------------------------------------------------------------------
    % S4 index calculation
    
%     u_noise_int = spline(slta./1e3, u_noise, slta_int);
%     
%     u_signal_int = spline(slta./1e3, u_signal, slta_int);
    
%     ix = slta >= 40e3 & slta <= 580e3; % modified 18/07/23

    ix = slta >= 40e3;
    
    sci.snr = snr_interpolated(ix);
        
    sci.slta = slta(ix);
    
    sci.u_sq = uZverevBP(ix).*conj(uZverevBP(ix));
    
    sci.sigma_0 = sigma_0(ix);
    
    sci.fsWOP = fsWOP;
    
    sci.fs = fsMeas;
    
    sci.dirSave = data.dirSave;
    
    sci.fileName = strrep(figName, 'noise', 'S4index');
    
    sci.date = data.date;
    
%     [sci] = s4index_w_noise(sci); % uncomment to compute S4 index
    
% %     return;
    
    % ---------------------------------------------------------------------

    figure(171); plot(abs(uZverevBP), slta./1e3, 'b', 'DisplayName', 'u_{total}'); hold on
    
    plot(abs(u_signal), slta./1e3, 'r', 'DisplayName', 'u_{signal}');
    
    ylabel('SLTA (km)');
    
    xlabel('abs(u)');
    
    legend show;
    
    set(gca, 'FontName', 'Courier');
    
    ylim([40 600])
    
    yticks([100:100:600]);
    
    xlim([0.4 1.4])
    
%     t = datetime('now');
%     
%     t.Format = 'yyyy-MM-dd''T''HH_mm_ss';
    
    set(gcf, 'Position', [445.4000  426.2000  314.8000  420.0000]);
    
    fullName = strrep(fullfile(figFolder, figName), '.mat', ['_proc' data.date '.png']);
    
    save(fullName, 'slta_int', 'sigma_0_int');
    
    figName = strrep(figName, '.mat', ['_proc' data.date '.png']);
    
    %saveas(gcf, fullfile(figFolder, figName), 'png'); %uncomment to save
    %noise figure
    
    close all;
   
end