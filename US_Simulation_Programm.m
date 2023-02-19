%Material Simulation - Thurner Stephan Jannuary 2023
clear all
close all
%eps0 = 8.854*10^-12 F/m
%enter diamter of sample (assumption round transdcuer, can also be changed)
d = 0.00367;
%______add matrial that can be used for simulation______
PZT = struct('d',0.015,'t',150*10^-6,'rho',7800,'v',4674*(1+1i*0.012),'h',2.19*(1+1i*0.029)*10^9,'eps',1.06*(1-1i*0.053)*10^-8);
PVDFTrFE = struct('d',d,'t',28*10^-6,'rho',1880,'v',2400*(1),'h',-4.7*(1)*10^9,'eps',4.692*(1-1i*0.14)*10^-11);
PVDF = struct('d',0.010,'t',28*10^-6,'rho',1780,'v',2260*(1),'h',-2.6*(1)*10^9,'eps',5.489*(1-1i*0.25)*10^-11);
water = struct('d',d,'t',1*10^-6,'rho',1000,'v',1500);
air = struct('d',0.015,'t',0.001,'rho',1.24,'v',344);
silverInk = struct('d',d,'t',0.5*10^-6,'rho',3400,'v',3600);
silverepoxy = struct('d',0.010,'t',6*10^-6,'rho',2710,'v',1900);
polyimide = struct('d',d,'t',50*10^-6,'rho',1421,'v',2414);
mylar = struct('d',d,'t',13*10^-6,'rho',1180,'v',2540);
epoxy = struct('d',d,'t',28*10^-6,'rho',1100,'v',2200*(1+1i*0.05));
steel = struct('d',0.001,'t',150*10^-6,'rho',7890,'c',2.645*(1+1i*0.002)*10^9,'v',5790*(1+1i*0.001));
Ni = struct('d',d,'t',0.15*10^-6,'rho',8902,'v',6040);
Cu = struct('d',d,'t',0.4*10^-6,'rho',8960,'v',4760);
silicon = struct('d',d,'t',12*10^-6,'rho',1470,'v',960);
%______left and right accoustic port
ZL = 0;
ZR = 0;
% ZR = (water.rho)*(water.v)*(((water.d)^2)*pi/4)

%______add posible material setups that will be used for simulation_____
%US = struct('piezo',PVDF);
%US = struct('backing',steel,'piezo',PZT);
%US = struct('bottomEL',silverInk,'piezo',PVDFTrFE,'topEL',silverInk,'backing',polyimide,'water',water);
US = struct('Ni_b',Ni,'Cu_b',Cu,'piezo',PVDFTrFE,'Cu_f',Cu,'Ni_f',Ni);
% US = struct('bottomEL',silverInk,'piezo',PVDFTrFE,'topEL',silverInk,'mylar',polyimide);
%US = struct('piezo',PZT);
% US = struct('piezo',PVDFTrFE,'CuBacking',Cu);
% US = struct('piezo',PVDFTrFE);
% US = struct('piezo',PVDFTrFE,'Adhnesive',epoxy,'SiliconBacking',silicon);

%______do not change, ==> matrix for the caluclation is created
name = fieldnames(US).';
US_shape = size(fieldnames(US));
%matrix for data saving
Zi = zeros(3,US_shape(1,1)+4);
%add terminating resistor
Zi(1,1) = ZL;
Zi(1,US_shape(1,1)+4) = ZR;
%get postition of piezo element
piezo_pos = find(strcmp(name,'piezo'));
%get piezo properties
C0 = (US.(string(name(piezo_pos))).eps)*(((US.(string(name(piezo_pos))).d)^2)*pi/4)/(US.(string(name(piezo_pos))).t);
N = C0*US.(string(name(piezo_pos))).h;
%calculate impedance Z0 of each medium >Z0=rho*A*v<
for i = 1:(US_shape(1,1)+2)
    if i<=piezo_pos
        Zi(1,i+1) = (US.(string(name(i))).rho)*(US.(string(name(i))).v)*(((US.(string(name(i))).d)^2)*pi/4);
    end
    if i>piezo_pos+1
        Zi(1,i+1) = (US.(string(name(i-2))).rho)*(US.(string(name(i-2))).v)*(((US.(string(name(i-2))).d)^2)*pi/4);
    end
end

%______uncomment functions to create desired plots______ 
%>>>>(f_min,f_max,f_step) in MHz
%____shows a graph that describes the impedance over the frequency
% plot_freq_domain(1,240,0.1,Zi,US,US_shape,name,piezo_pos,C0,N)
%____analyse each parameter of simulation
%  plot_freq_sub(0,0,80,0.1,Zi,US,US_shape,name,piezo_pos,C0,N)
%>>>(f_mean,f_min,f_max,f_step) in MHz and if f_mean == 0 > plot resonance freq<<<
%____analyse multiple VNA analyed samples%plot_freq_sub_VNA(1,200,0.01,Zi,US,US_shape,name,piezo_pos,C0,N)
% plot_freq_sub_VNA(1,140,0.1,Zi,US,US_shape,name,piezo_pos,C0,N)
%____analyze the VNA samples each plot separate
% plot_freq_sep_VNA(1,160,0.1,Zi,US,US_shape,name,piezo_pos,C0,N)
%____find a suitable filter structure for a VNA analyzed sample >>pure piezo<<
% plot_matching(1,60,0.01,Zi,US,US_shape,name,piezo_pos,C0,N)
%____compare the real filter to the theoretic setup >>VNA pure piezo, matched<<
% plot_matching_compare("LL",15,55,0.01,Zi,US,US_shape,name,piezo_pos,C0,N)

%--------------------------------------------------function-------------------------------------------------
%plot_matching_in_comparisation_to_realtiy
function plot_matching_compare(filter,f_min,f_max,f_step,Zi,US,US_shape,name,piezo_pos,C0,N)
    VNA = ["C:\export\MA_Verif\PL2.csv";
           "C:\export\MA_Verif\ML2.csv"];
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot frequency domin%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = f_min*10^6:f_step*10^6:f_max*10^6;
    x_plot = x./10^6;
    Zcal = zeros(size(x));    
    for i = 1:length(x)
        impedance = calculate_Z(x(i),Zi,US,US_shape,name,piezo_pos,C0,N);
        Zcal(i) = impedance(1,piezo_pos+2);
        gamma_gen(i) = (Zcal(i)-50)/(Zcal(i)+50);
        VSWR(i) = (1+abs(gamma_gen(i)))/(1-abs(gamma_gen(i)));
    end
    %set up theoretic adaption to reality
    Zcalsmooth = smooth_theo_Z(Zcal,0.02);
    Zabs = abs(Zcal);
    %get minimums of array
    [pks,locs] = findpeaks(-Zabs);
    min_index = x(locs)./10^6;
    %get name of assabeld transducer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%matching%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xo, RLo, PH1o, RSo, SWRo, XSo, Zo] = get_VNA_data(string(VNA(1)));
    [xm, RLm, PH1m, RSm, SWRm, XSm, Zm] = get_VNA_data(string(VNA(2)));
    %only take the first value
    f0 = min_index;
    if length(f0) > 1
        f0 = f0(1);
    end
    %find matching position of the VNA
    res_index=0;
    x_mega = xo./10^6;
    for i = 1:length(xo)
        if x_mega(i) > f0
            res_index = i;
            break;
        end
    end
    %add data to match, there u can choose the dataset that should be
    %matched either simulation or VNA real data
    xmatch = xo;
    Immatch = XSo;
    Rematch = RSo;
    f0_VNA = xmatch(res_index);
    Im_f0 = Immatch(res_index);
    Re_f0 = Rematch(res_index);
    %create plot matrix inculding all important data
    PlotMatch = ones(5,length(xmatch));
    for i=1:length(xmatch)
        PlotMatch(1,i) = xmatch(i);
        PlotMatch(2,i) = complex(Rematch(i),Immatch(i));
        PlotMatch(4,i) = (PlotMatch(2,i)-50)/(PlotMatch(2,i)+50);
    end
    PlotMatchReal = ones(3,length(xm));
    for i=1:length(xm)
        PlotMatchReal(1,i) = xm(i);
        PlotMatchReal(2,i) = complex(RSm(i),XSm(i));
        PlotMatchReal(3,i) = (PlotMatchReal(2,i)-50)/(PlotMatchReal(2,i)+50);
    end
    %calcualte paramters for choosen network
    EM = match_Z(Im_f0,Re_f0,f0_VNA,filter);
    disp('Matching network: '+string(filter));
    display_matching(EM);
    %sweep over frequency range
    for n=1:length(xmatch)
        %sweep each element
        %>>>sweep_matching(Im_f0,Re_f0,f0_VNA,EM)<<<;
        PlotMatch(3,n) = sweep_matching(imag(PlotMatch(2,n)),real(PlotMatch(2,n)),PlotMatch(1,n),EM);
    end
    for m=1:length(xmatch)
        %sweep each element
        %calucalte return loss
        PlotMatch(5,m) = (PlotMatch(3,m)-50)/(PlotMatch(3,m)+50);
    end
    %%%%%%%%%%%%%%%%plot data%%%%%%%%%%%%%%%
    figure('Name','Matching Network: Compare to reality','WindowState','maximized');
    subplot(2,2,1)
    hold on
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [0 500], 'Color',[1 0 1],'LineStyle','-.','LineWidth',2,'DisplayName',string(round(min_index(p),0))+"MHz (sim)");
    end
    %plot(x_plot,abs(Zcalsmooth),'DisplayName','Z theo');
    plot(PlotMatch(1,:)./10^6,abs(PlotMatch(2,:)),'Color',[0 0 1],'LineWidth',2,'DisplayName','VNA pure piezo');
    plot(PlotMatchReal(1,:)./10^6,abs(PlotMatchReal(2,:)),'Color',[1 0.2 0],'LineWidth',2,'DisplayName','VNA matched');
    plot(PlotMatch(1,:)./10^6,abs(PlotMatch(3,:)),'Color',[0 1 0],'LineStyle','-','LineWidth',2,'DisplayName',string(filter+' simulation'));
    title("Impedance");
    xlabel('Frequency [MHz]');
    ylabel('Impedance Z [Ohm]');
    legend show
    xlim([10,f_max]);
    %set(gca, 'YScale', 'log')
    grid on
    ylim([0 400])
    xlim([15 55])
    hold off
    %RL dB
    subplot(2,2,4)
    hold on
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [-80 0], 'Color',[1 0 1],'LineStyle','-.','LineWidth',2,'DisplayName',string(round(min_index(p),0))+"MHz (sim)");
    end
    plot(PlotMatch(1,:)./10^6,10*log(abs(PlotMatch(4,:))),'Color',[0 0 1],'LineWidth',2,'DisplayName','VNA pure piezo');
    plot(PlotMatchReal(1,:)./10^6,10*log(abs(PlotMatchReal(3,:))),'Color',[1 0.2 0],'LineWidth',2,'DisplayName','VNA matched');
    plot(PlotMatch(1,:)./10^6,10*log(abs(PlotMatch(5,:))),'Color',[0 1 0],'LineStyle','-','LineWidth',2,'DisplayName',string(filter+' simulation'));
    title("Returnloss dB");
    xlabel('Frequency [MHz]');
    ylabel('Returnloss RL [dB]');
    legend show
    %set(gca, 'YScale', 'log')
    grid on
    ylim([-40 0])
    xlim([15 55])
    hold off
    %RS&XS
    subplot(2,2,2)
    hold on
    yyaxis left
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [0 500], 'Color',[1 0 1],'LineStyle','-.','LineWidth',2,'DisplayName',string(round(min_index(p),0))+"MHz (sim)");
    end
    plot(PlotMatch(1,:)./10^6,real(PlotMatch(2,:)),'Color',[0 0 1],'LineStyle','-','LineWidth',2,'DisplayName','VNA pure piezo [RS]');
    plot(PlotMatchReal(1,:)./10^6,real(PlotMatchReal(2,:)),'Color',[1 0.2 0],'LineStyle','-','LineWidth',2,'DisplayName','VNA matched [RS]');
    plot(PlotMatch(1,:)./10^6,real(PlotMatch(3,:)),'Color',[0 1 0],'LineStyle','-','LineWidth',2,'DisplayName',string(filter+' simulation [RS]'));
    yyaxis right
    plot(PlotMatch(1,:)./10^6,imag(PlotMatch(2,:)),'Color',[0 0 1],'LineStyle','-.','LineWidth',2,'DisplayName','VNA pure piezo [XS]');
    plot(PlotMatchReal(1,:)./10^6,imag(PlotMatchReal(2,:)),'Color',[1 0.2 0],'LineStyle','-.','LineWidth',2,'DisplayName','VNA matched [XS]');
    plot(PlotMatch(1,:)./10^6,imag(PlotMatch(3,:)),'Color',[0 1 0],'LineStyle','-.','LineWidth',2,'DisplayName',string(filter+' simulation [XS]'));
    yyaxis left
    title("RS & XS impedance");
    xlabel('Frequency [MHz]');
    ylabel('real impdeance RS [Ohm]');
    legend show
    grid on
    ylim([0 200])
    yyaxis right
    ylabel('blind impdeance XS [Ohm]');
    legend show
    grid on
    xlim([15,55]);
    %set(gca, 'YScale', 'log')
    hold off
end
%--------------------------------------------------function-------------------------------------------------
%plot_matching_suitable_filter_structure
function plot_matching(f_min,f_max,f_step,Zi,US,US_shape,name,piezo_pos,C0,N)
    VNA = ["C:\export\MA_Verif\PS4.csv";
           "C:\export\MA_Verif\PM.csv";
           "C:\export\MA_Verif\PL2.csv";];
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot frequency domin%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = f_min*10^6:f_step*10^6:f_max*10^6;
    x_plot = x./10^6;
    yz = zeros(size(x));
    for i = 1:length(x)
        impedance = calculate_Z(x(i),Zi,US,US_shape,name,piezo_pos,C0,N);
        yz(i) = impedance(1,piezo_pos+2);
        gamma_gen(i) = (yz(i)-50)/(yz(i)+50);
        VSWR(i) = (1+abs(gamma_gen(i)))/(1-abs(gamma_gen(i)));
    end
    y = abs(yz);
    %get minimums of array
    [pks,locs] = findpeaks(-y);
    min_index = x(locs)./10^6;
    %get name of assabeld transducer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%matching%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %select data set
    [xo, RLo, PH1o, RSo, SWRo, XSo, Zo] = get_VNA_data(string(VNA(3)));
    %only take the first value
    f0 = min_index;
    if length(f0) > 1
        f0 = f0(1);
    end
    %find matching position of the VNA
    res_index=0;
    x_mega = xo./10^6;
    for i = 1:length(xo)
        if x_mega(i) > f0
            res_index = i;
            break;
        end
    end
    %add data to match, there u can choose the dataset that should be
    %matched either simulation or VNA real data
    xmatch = xo;
    Immatch = XSo;
    Rematch = RSo;
    f0_VNA = xmatch(res_index);
    Im_f0 = Immatch(res_index);
    Re_f0 = Rematch(res_index);
    %find suitable matching network
    MatchNet = ["resS";"resP";"LS";"LP";"LL";"LLS"];
    %create plot matrix inculding all important data
    PlotMatch = ones(17,length(xmatch));
    for i=1:length(xmatch)
        PlotMatch(1,i) = xmatch(i);
        %>>>reality or stick to a value<<<
        PlotMatch(2,i) = complex(Rematch(i),Immatch(i));
        %PlotMatch(2,i) = complex(Re_f0,Im_f0);
        PlotMatch(10,i) = (PlotMatch(2,i)-50)/(PlotMatch(2,i)+50);
    end
    for i=1:6
        %calcualte paramters for each network
        EM = match_Z(Im_f0,Re_f0,f0_VNA,MatchNet(i));
        disp('Matching network: '+string(MatchNet(i)));
        display_matching(EM);
        %sweep over frequency range
        for n=1:length(xmatch)
            %sweep each element
            %>>>sweep_matching(Im_f0,Re_f0,f0_VNA,EM)<<<;
            PlotMatch(2+i,n) = sweep_matching(imag(PlotMatch(2,n)),real(PlotMatch(2,n)),PlotMatch(1,n),EM);
        end
        for m=1:length(xmatch)
            %sweep each element
            %calucalte return loss
            PlotMatch(10+i,m) = (PlotMatch(2+i,m)-50)/(PlotMatch(2+i,m)+50);
        end
    end
    %%%%%%%%%%%%%%%%plot data%%%%%%%%%%%%%%%
    figure('Name','Matching Network: Suitable filter structures ','WindowState','maximized');
    subplot(2,2,1)
    hold on
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [0 800], 'Color',[1 0 1],'LineStyle','-.','LineWidth',2,'DisplayName',string(round(min_index(p),0))+"MHz (sim)");
    end
    plot(PlotMatch(1,:)./10^6,abs(PlotMatch(2,:)),'LineWidth',2,'DisplayName','VNA real');
    for i=1:6
        plot(PlotMatch(1,:)./10^6,abs(PlotMatch(2+i,:)),":",'LineWidth',2,'DisplayName',(MatchNet(i)));
    end
    title("Impedance");
    xlabel('Frequency [MHz]');
    ylabel('Impedance Z [Ohm]');
    legend show
    xlim([f_min,f_max]);
    %set(gca, 'YScale', 'log')
    grid on
    ylim([0 640])
    xlim([15 55])
    hold off
    %RL dB
    subplot(2,2,4)
    hold on
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [-80 0], 'Color',[1 0 1],'LineStyle','-.','LineWidth',2,'DisplayName',string(round(min_index(p),0))+"MHz (sim)");
    end
    plot(PlotMatch(1,:)./10^6,10*log(abs(PlotMatch(10,:))),'LineWidth',2,'DisplayName','VNA real');
    for i=1:6
        plot(PlotMatch(1,:)./10^6,10*log(abs(PlotMatch(10+i,:))),":",'LineWidth',2,'DisplayName',(MatchNet(i)));
    end
    title("Returnloss dB");
    xlabel('Frequency [MHz]');
    ylabel('Returnloss RL [dB]');
    legend show
    xlim([f_min,f_max]);
    %set(gca, 'YScale', 'log')
    grid on
    ylim([-40 0])
    xlim([15 55])
    hold off
    %RS&XS
    subplot(2,2,2)
    hold on
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [0 800], 'Color',[1 0 1],'LineStyle','-.','LineWidth',2,'DisplayName',string(round(min_index(p),0))+"MHz (sim)");
    end
    plot(PlotMatch(1,:)./10^6,real(PlotMatch(2,:)),'LineWidth',2,'DisplayName','VNA real');
    for i=1:6
        plot(PlotMatch(1,:)./10^6,real(PlotMatch(2+i,:)),":",'LineWidth',2,'DisplayName',(MatchNet(i)));
    end
    title("Real impedance");
    xlabel('Frequency [MHz]');
    ylabel('real impdeance RS [Ohm]');
    legend show
    xlim([f_min,f_max]);
    %set(gca, 'YScale', 'log')
    grid on
    %ylim([-30 0])
    xlim([15 55])
    hold off
end
%--------------------------------------------------function-------------------------------------------------
%plot_multiple_data_of_the_VNA
function plot_freq_sep_VNA(f_min,f_max,f_step,Zi,US,US_shape,name,piezo_pos,C0,N)
    VNA = ["C:\export\Euromat\20-200\T61.csv";
           "C:\export\Euromat\20-200\T62.csv";
           "C:\export\Euromat\20-200\T63.csv";
           "C:\export\Euromat\20-200\T64.csv";
           "C:\export\Euromat\20-200\T66.csv";
           "C:\export\Euromat\20-200\T67.csv"];
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot frequency domin%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = f_min*10^6:f_step*10^6:f_max*10^6;
    x_plot = x./10^6;
    yz = zeros(size(x));
    for i = 1:length(x)
        impedance = calculate_Z(x(i),Zi,US,US_shape,name,piezo_pos,C0,N);
        yz(i) = impedance(1,piezo_pos+2);
        gamma_gen(i) = (yz(i)-50)/(yz(i)+50);
        VSWR(i) = (1+abs(gamma_gen(i)))/(1-abs(gamma_gen(i)));
    end
    y = abs(yz);
    %get minimums of array
    [pks,locs] = findpeaks(-y);
    [pks1,locs1] = findpeaks(-y);
    [pks2,locs2] = findpeaks(y);
    min_index1 = x(locs1)./10^6;
    min_index2 = x(locs2)./10^6;
    min_index = (min_index2-min_index1)/2+min_index1;
    plotr = min_index(1);
    f_mean = min_index;
    val_min = min(y);
    val_max = max(y);
    %get name of assabeld transducer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Z%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Name','Impedance Z of the VNA samples','WindowState','maximized');
    hold on
    plot(x_plot,y,'b:','DisplayName','Impedance Z Simulation')
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [6 val_max], 'Color',[1 0 1],'LineStyle','--','DisplayName',string(round(min_index(p),0))+" MHz ("+string(round(-pks(p)))+" Ohm)");
    end
    for i=1:length(VNA)
        [xo, RLo, PH1o, RSo, SWRo, XSo, Zo] = get_VNA_data(string(VNA(i)));
        plot(xo./10^6,Zo);
    end
    title("Impedance of VNA");
    xlabel('Frequency [MHz]');
    ylabel('Impedance Z [Ohm]');
    legend show
    set(gca, 'YScale', 'log')
    grid on
    ylim([0 400])
%     xlim([5 60])
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Name','Real resistance RS and Blind resistance XS of VNA samples','WindowState','maximized');
    hold on
    for i=1:length(VNA)
        [xo, RLo, PH1o, RSo, SWRo, XSo, Zo] = get_VNA_data(string(VNA(i)));
        yyaxis left
        plot(xo./10^6,RSo);
        yyaxis right
        plot(xo./10^6,XSo);
    end   
    yyaxis left
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [0 400], 'Color',[1 0 1],'LineStyle','--','DisplayName',string(round(min_index(p),0))+" MHz");
    end
    hold off
    title("Real resistance of VNA");
    xlabel('Frequency [MHz]');
    ylabel('Impedance RS [Ohm]');
    ylim([0 400])
    yyaxis right
    ylabel('Impedance XS [Ohm]');
    legend show
    %set(gca, 'YScale', 'log')
%     xlim([5 60])
    grid on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Name','Phase of VNA samples','WindowState','maximized');
    hold on
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [-2000 val_max], 'Color',[1 0 1],'LineStyle','--','DisplayName',string(min_index(p))+" MHz ("+string(round(-pks(p)))+" Ohm)");
    end
    for i=1:length(VNA)
        [xo, RLo, PH1o, RSo, SWRo, XSo, Zo] = get_VNA_data(string(VNA(i)));
        plot(xo./10^6,PH1o);
    end
    title("Phase of VNA");
    xlabel('Frequency [MHz]');
    ylabel('Phase');
    legend show
    %set(gca, 'YScale', 'log')
    grid on
    ylim([-200 200])
%     xlim([5 80])
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Name','Return loss of VNA samples','WindowState','maximized');
    hold on
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [-2000 val_max], 'Color',[1 0 1],'LineStyle','--','DisplayName',string(round(min_index(p),0))+" MHz ("+string(round(-pks(p)))+" Ohm)");
    end
    for i=1:length(VNA)
        [xo, RLo, PH1o, RSo, SWRo, XSo, Zo] = get_VNA_data(string(VNA(i)));
        plot(xo./10^6,RLo);
    end
    title("Returnloss of VNA in dB");
    xlabel('Frequency [MHz]');
    ylabel('Returnloss [dB]');
    legend show
    %set(gca, 'YScale', 'log')
    grid on
    ylim([-12 0])
%     xlim([5 80])
    hold off
end
%--------------------------------------------------function-------------------------------------------------
%plot data combined with the VNA
function plot_freq_sub_VNA(f_min,f_max,f_step,Zi,US,US_shape,name,piezo_pos,C0,N)
  VNA = ["C:\export\standard\s_200.csv"];
    
    %%% get data from VNA
    [x_VNA, RL_VNA, PH_VNA, RS_VNA, SWR_VNA, XS_VNA, Z_VNA] = get_VNA_data(VNA(1));
    x_VNA_plot = x_VNA./10^6;
    % get data from VNA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot frequency domin%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = f_min*10^6:f_step*10^6:f_max*10^6;
    x_plot = x./10^6;
    yz = zeros(size(x));
    for i = 1:length(x)
        impedance = calculate_Z(x(i),Zi,US,US_shape,name,piezo_pos,C0,N);
        yz(i) = impedance(1,piezo_pos+2);
        gamma_gen(i) = (yz(i)-50)/(yz(i)+50);
        VSWR(i) = (1+abs(gamma_gen(i)))/(1-abs(gamma_gen(i)));
    end
    y = abs(yz);
    %get minimums of array
    [pks,locs] = findpeaks(-y);
    [pks1,locs1] = findpeaks(-y);
    [pks2,locs2] = findpeaks(y);
    min_index1 = x(locs1)./10^6;
    min_index2 = x(locs2)./10^6;
    min_index = (min_index2-min_index1)/2+min_index1;
    plotr = min_index(1);
    f_mean = min_index;
    val_min = min(y);
    val_max = max(y);
    %get name of assabeld transducer
    leg_name="Freq domain";
    %plot frequency behaviour
    figure('Name','Freq domain and wave propagation','WindowState','maximized');
    subplot(2,2,1)
    hold on;
    plot(x_plot,y,'DisplayName','Impedance Z Simulation')
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [1 4000], 'Color',[1 0 1],'LineStyle','-','DisplayName',string(round(min_index(p),0))+" MHz ("+string(round(-pks(p)))+" Ohm)");
    end
    %plot VNA data
    for i=1:length(VNA)
        [xo, RLo, PH1o, RSo, SWRo, XSo, Zo] = get_VNA_data(VNA(i));
        plot(xo./10^6,Zo,'LineStyle','--');
    end
    title(leg_name);
    xlabel('Frequency [MHz]');
    ylabel('Impedance Z [Ohm]');
    legend show
    set(gca, 'YScale', 'log')
    grid on
    ylim([10 4000])
    xlim([1 140])
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%impedance%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,2)
    %get name of assabeld transducer
    leg_name="RS & XS of VNA samples";
    hold on
    for i=1:length(VNA)
        [xo, RLo, PH1o, RSo, SWRo, XSo, Zo] = get_VNA_data(VNA(i));
        yyaxis left
        plot(xo./10^6,RSo,'DisplayName',string('RS'+string(i)));
        yyaxis right
        plot(xo./10^6,XSo,'DisplayName',string('XS'+string(i)));
    end
    yyaxis left
    for p = 1:1
        line([min_index(p) min_index(p)], [0 1000], 'Color',[1 0 1],'LineStyle','--','DisplayName',string(round(min_index(p),0))+" MHz");
    end
    title(leg_name);
    xlabel('Frequency [MHz]');
    yyaxis left
    ylabel('Impedance RS [Ohm]');
    yyaxis right
    ylabel('Impedance XS [Ohm]');
    legend show
    xlim([1,140])
    grid on
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%VNA measurements%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,1,2)
    %get name of assabeld transducer
    leg_name="VNA measurements normalized";
    hold on
    
    VNAE = 400;
    for i=1:1%length(VNA)
        [xo, RLo, PH1o, RSo, SWRo, XSo, Zo] = get_VNA_data(VNA(i));
        if VNAE == 0
            plot(xo./10^6,normalize(Zo,'range'),'DisplayName','Z');
            plot(xo./10^6,normalize(RSo,'range'),'DisplayName','RS');
            plot(xo./10^6,normalize(XSo,'range'),'DisplayName','XS');
            plot(xo./10^6,normalize(PH1o,'range'),'DisplayName','PH');
            plot(xo./10^6,normalize(RLo,'range'),'DisplayName','RL');
            plot(xo./10^6,normalize(SWRo,'range'),'DisplayName','SWR');
        else
            plot(xo(1:VNAE)./10^6,normalize(Zo(1:VNAE),'range'),'DisplayName','Z');
            plot(xo(1:VNAE)./10^6,normalize(RSo(1:VNAE),'range'),'DisplayName','RS');
            plot(xo(1:VNAE)./10^6,normalize(XSo(1:VNAE),'range'),'DisplayName','XS');
            plot(xo(1:VNAE)./10^6,normalize(PH1o(1:VNAE),'range'),'DisplayName','PH');
            plot(xo(1:VNAE)./10^6,normalize(RLo(1:VNAE),'range'),'DisplayName','RL');
            plot(xo(1:VNAE)./10^6,normalize(SWRo(1:VNAE),'range'),'DisplayName','SWR');
        end
    end
    for p = 1:length(min_index)
        line([min_index(p) min_index(p)], [0 1], 'Color',[1 0 1],'LineStyle','--','DisplayName',"Simulation");
    end
    hold off
    title(leg_name);
    xlabel('Frequency [MHz]');
    ylabel('VNA measurements [0-1]');
    legend show
    grid on
    xlim([1 80]);
end
%--------------------------------------------------function-------------------------------------------------
%plot multiple data
function plot_freq_sub(f_mean,f_min,f_max,f_step,Zi,US,US_shape,name,piezo_pos,C0,N)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot frequency domin%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_sym = 0; %enable the symetry of the wave propagation in the piezo
    piezo_in_middle = 0; %choose if the wave is centert at the piezo
    x = f_min*10^6:f_step*10^6:f_max*10^6;
    x_plot = x./10^6;
    y = zeros(size(x));
    for i = 1:length(x)
        impedance = calculate_Z(x(i),Zi,US,US_shape,name,piezo_pos,C0,N);
        y(i) = abs(impedance(1,piezo_pos+2));
    end
    %get minimums of array
    [pks,locs] = findpeaks(-y);
    [pks1,locs1] = findpeaks(-y);
    [pks2,locs2] = findpeaks(y);
    min_index1 = x(locs1)./10^6;
    min_index2 = x(locs2)./10^6;
    min_index = (min_index2-min_index1)/2+min_index1;
    plotr = min_index(1);
    f_mean = min_index;
%     [pks,locs] = findpeaks(-y);
%     min_index = x(locs)./10^6;
    val_min = min(y);
    val_max = max(y);
    %get name of assabeld transducer
    leg_name="Freq domain: ";
    for n = 1:US_shape(1,1)
        leg_name = leg_name + string(name(n))+ "("+string((US.(string(name(n))).t)*10^6)+"um) ";
    end
    %plot frequency behaviour
    figure('Name','Freq domain and wave propagation','WindowState','maximized');
    subplot(2,2,1)
    plot(x_plot,y,'DisplayName','Impedance Z')
    hold on;
    for p = 1:length(min_index)
        min_index(p);
        line([min_index(p) min_index(p)], [2 val_max], 'Color',[1 0 1],'LineStyle','--','DisplayName',string(round(min_index(p),0))+" MHz");
    end
    title(leg_name);
    xlabel('Frequency [MHz]');
    ylabel('Impedance Z [Ohm]');
    legend show
    set(gca, 'YScale', 'log')
    grid on
    ylim([10 4*10^3])
    xlim([0 80])
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot geometry%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot at resoncance frequences
    if f_mean == 0
        f_mean = x_plot(locs);
    end
    subplot(2,1,2)
    %get name of assabeld transducer
    leg_name="Z0 | RL: ";
    dist_name = [];
    for n = 1:US_shape(1,1)
        leg_name = leg_name + string(name(n))+ "("+string((US.(string(name(n))).t)*10^6)+"um) ";
        dist_name = [dist_name, US.(string(name(n))).t];
    end
    dist_name = [0.1*sum(dist_name) dist_name 0.1*sum(dist_name)];
    %create distance matrix EN-Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %>>>enter resolution<<<
    y_en = zeros(length(dist_name),1000);
    x = linspace(0,1,length(y_en))*sum(dist_name);
    el_inx = 1;
    dist_inc = 0;
    for i = 1:length(x)
        if x(i) < dist_name(el_inx)+dist_inc
            y_en(el_inx,i) = 1;
        elseif el_inx < length(dist_name)
            dist_inc = dist_inc+dist_name(el_inx);
            el_inx = el_inx+1;
            y_en(el_inx,i) = 1;
        end
    end
    %get impdance data
    y_Z0_dist = [];
    y_Zi_dist = [];
    for res_i = 1:length(f_mean)
        impedance = calculate_Z(f_mean(res_i)*10^6,Zi,US,US_shape,name,piezo_pos,C0,N);
        y_Z0 = abs(impedance(1,:));
        y_Zi = abs(impedance(2,:));
        %remove double dataset
        y_Z0(piezo_pos+1) = [];
        y_Z0(piezo_pos+1) = [];
        y_Zi(piezo_pos+1) = [];
        y_Zi(piezo_pos+2) = [];
        %create impedance vector
        y_Z0_dist = [y_Z0_dist; (y_en'*y_Z0')'];
        y_Zi_dist = [y_Zi_dist; (y_en'*y_Zi')'];
    end
    %find min and max of matrix
    val_min = min(y_Z0_dist(:));
    val_max = max(y_Z0_dist(:));
    if val_min < 10^-2
        val_min = 10^-2;
    end
    %calucalte markers for distance
    material_dist = diff((y_en'*(1:1:length(dist_name))')');
    material_inx = find(material_dist==1);
    %plot data
    hold on;
    yyaxis left
    plot(x*10^6,y_Z0_dist(1,:),'b','DisplayName','Impedance Z0');
%     for res_i = 1:length(f_mean)
%         plot(x*10^6,y_Zi_dist(res_i,:),'b--','DisplayName','Impedance Zin: '+string(round(f_mean(res_i),0))+'MHz');
%     end
    for i = 1:length(material_inx)
        line([x(material_inx(i))*10^6 x(material_inx(i))*10^6], [val_min val_max*1.1],'Color',[0 1 0],'LineStyle','-.','DisplayName',string(round(x(material_inx(i)).*10^6))+" um");
    end
    hold off;
    title(leg_name);
    xlabel('Distance [um]');
    yyaxis left
    ylabel('Impedance Z [Ohm]');
    legend show
    ylim([1,val_max*1.1]);
%     set(gca, 'YScale', 'log')
    grid on;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%return coefficient%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get impedance data
    y_r0_dist = [];
    for res_i = 1:length(f_mean)
        impedance = calculate_Z(f_mean(res_i)*10^6,Zi,US,US_shape,name,piezo_pos,C0,N);
        %y_r0 = abs(impedance(3,:));
        y_r0 = impedance(3,:);
        y_r0(piezo_pos+1) = [];
        y_r0(piezo_pos+1) = [];
        %create reflection vector
        y_r0_array = (y_en'*y_r0')';
        y_r0_array(end) = y_r0_array(end-1);
        y_r0_dist = [y_r0_dist; y_r0_array];
    end
    yyaxis right
    hold on;
    for res_i = 1:length(f_mean)
        plot(x*10^6,y_r0_dist(res_i,:),'Color',[1 0 0],'LineStyle','-','DisplayName','Return coefficient: '+string(round(f_mean(res_i),0))+'MHz');
    end
    ylabel('Return coefficient');
    ylim([-1.2 0.2]);
    legend show
%     set(gca, 'YScale', 'log')
    grid on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wave propagation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,2)
    %get name of assabeld transducer
    %disable center of piezo
    if piezo_in_middle==0
        piezo_pos = 1;
    end
    leg_name="Wavelength: ";
    for n = 1:US_shape(1,1)
        leg_name = leg_name + string(name(n))+ "("+string((US.(string(name(n))).t)*10^6)+"um) ";
    end
    y_sin_dist = [];
    y_sin_sym_dist = [];
    %calucalte for each resonance freq
    for res_i = 1:length(f_mean)
        %create beta array to calculate the
        beta_array = [];
        phase_array = [];
        for p = 1:US_shape(1,1)
            beta_array = [beta_array abs(wave_const(f_mean(res_i)*10^6,US.(string(name(p))).v))];
            phase_array = [phase_array abs(wave_const(f_mean(res_i)*10^6,US.(string(name(p))).v))*US.(string(name(p))).t];
        end
        %phase_array = phase_array*1.5;
        %calcualte points in material
        material_array = [];
        for p = 1:length(dist_name)
            material_array = [material_array sum(y_en(p,:) == 1)];
        end
        %create sin for each
        y_sin = [];
        for p = 1:length(dist_name)
            if p == 1 || p == length(dist_name)
                y_sin = [y_sin sin(linspace(0,0,material_array(p)))];
            else
                y_sin = [y_sin sin(linspace(0,phase_array(p-1),material_array(p)))];
            end
        end
        y_sin = [y_sin y_sin(end)];
        %add to matrix for different frequencys
        y_sin_dist = [y_sin_dist; y_sin];
        %make a contious sinus with the center at the piezo
        %move the sin to the center of the piezo
        x_sin_sym = linspace(0,phase_array(piezo_pos),material_array(piezo_pos+1));
        if plot_sym==1
            for i = 1:length(x_sin_sym)
                val = sin(x_sin_sym);
                %move elemts from the front into the back to make it symetric
                if abs(val(1)) < abs(val(end))
                    x_sin_sym = [x_sin_sym (x_sin_sym(end)+(abs(x_sin_sym(2))-abs(x_sin_sym(1))))];
                    x_sin_sym(1) = [];
                else
                    break
                end
            end
        end
        %add x values to the end of the piezo but they are symetric
        for p = piezo_pos+2:length(dist_name)
            if p ~= length(dist_name)
                x_sin_sym = [x_sin_sym linspace(x_sin_sym(end), x_sin_sym(end)+phase_array(p-1),material_array(p))];
            else
                x_sin_sym = [x_sin_sym sin(linspace(0,0,material_array(p)))];
            end
        end
        %add x values to the front of the piezo
        for r = 1:piezo_pos
            %reverse order start from the center
            p = piezo_pos+1-r;
            if p ~= 1
                x_sin_sym = [linspace(x_sin_sym(1)-phase_array(p-1), x_sin_sym(1),material_array(p)) x_sin_sym];
            else
                x_sin_sym = [sin(linspace(0,0,material_array(p))) x_sin_sym];
            end
        end
        %convert to sin and correct length of array
        y_sin_sym = sin(x_sin_sym);
        y_sin_sym = [y_sin_sym y_sin_sym(end)];
        %add to matrix for different frequencys
        y_sin_sym_dist = [y_sin_sym_dist; y_sin_sym];
    end
    %plot sin
    hold on;
    val_min = -1;
    val_max = 1;
    for res_i = 1:length(f_mean)
%         plot(x*10^6,y_sin_dist(res_i,:),'--','DisplayName','Invd.');
        plot(x*10^6,y_sin_sym_dist(res_i,:),'DisplayName','Sym.: '+string(round(f_mean(res_i),0))+'MHz');
    end
    for i = 1:length(material_inx)
        line([x(material_inx(i))*10^6 x(material_inx(i))*10^6], [val_min val_max], 'Color',[0 1 0],'LineStyle','-.','DisplayName',string(round(x(material_inx(i)).*10^6))+" um");
    end
    title(leg_name);
    xlabel('Distance [um]');
    ylabel('Propagation wavelength');
    legend show
    %set(gca, 'YScale', 'log')
    grid on
end
%--------------------------------------------------function-------------------------------------------------
%plot only impedance
function plot_freq_domain(f_min,f_max,f_step,Zi,US,US_shape,name,piezo_pos,C0,N)
    x = f_min*10^6:f_step*10^6:f_max*10^6;
    x_plot = x./10^6;
    y = zeros(size(x));
    for i = 1:length(x)
        impedance = calculate_Z(x(i),Zi,US,US_shape,name,piezo_pos,C0,N);
        y(i) = abs(impedance(1,piezo_pos+2));
    end
    %get minimums of array
    [pks,locs] = findpeaks(-y);
    [pks1,locs1] = findpeaks(-y);
    [pks2,locs2] = findpeaks(y);
    min_index1 = x(locs1)./10^6;
    min_index2 = x(locs2)./10^6;
    min_index = (min_index2-min_index1)/2+min_index1;
    plotr = min_index(1);
    f_mean = min_index;
%     [pks,locs] = findpeaks(-y);
%     min_index = x(locs)./10^6;
    val_min = min(y);
    val_max = max(y);
    %get name of assabeld transducer
    leg_name="Setup: ";
    for n = 1:US_shape(1,1)
        leg_name = leg_name + string(name(n))+ "("+string((US.(string(name(n))).t)*10^6)+"um) ";
    end
    %plot frequency behaviour
    figure('Name','Frequency domain of simulation','WindowState','maximized');
    hold on;
    plot(x_plot,y,'DisplayName','Impedance Z (sim)')
    for p = 1:length(min_index)
        min_index(p);
        line([min_index(p) min_index(p)], [val_min val_max], 'Color',[1 0 1],'LineStyle','--','DisplayName',string(min_index(p))+" MHz ("+string(round(-pks(p)))+" Ohm)");
    end
    hold off;
    title(leg_name);
    xlabel('Frequency [MHz]');
    ylabel('Impedance Z [Ohm]');
    ylim([0 10^3])
    set(gca, 'YScale', 'log')
    %set(gca, 'YScale', 'log')
    legend show
    grid on
end
%--------------------------------------------------function-------------------------------------------------
%get data from the VNA
function [freq, RL, PH, RS, SWR, XS, Z] = get_VNA_data(file_name)
    s=importdata(file_name);
    freq = (s(:,1))';
    RL = (s(:,2))';
    PH = (s(:,3))';
    RS = (s(:,4))';
    SWR = (s(:,5))';
    XS = (s(:,6))';
    Z = (s(:,7))';
end
%--------------------------------------------------function-------------------------------------------------
%display matching elements
function display_matching(EM)
    if EM(1,5) > 0
        disp('L_s0: '+string(EM(3,5))+' | '+string(EM(1,5)));
    elseif EM(1,5) < 0
        disp('C_s0: '+string(EM(3,5))+' | '+string(EM(1,5)));
    end
    if EM(2,4) > 0
        disp('L_p1: '+string(EM(4,4))+' | '+string(EM(2,4)));
    elseif EM(2,4) < 0
        disp('C_p1: '+string(EM(4,4))+' | '+string(EM(2,4)));
    end
    if EM(1,3) > 0
        disp('L_s1: '+string(EM(3,3))+' | '+string(EM(1,3)));
    elseif EM(1,3) < 0
        disp('C_s1: '+string(EM(3,3))+' | '+string(EM(1,3)));
    end
    if EM(2,2) > 0
        disp('L_p2: '+string(EM(4,2))+' | '+string(EM(2,2)));
    elseif EM(2,2) < 0
        disp('C_p2: '+string(EM(4,2))+' | '+string(EM(2,2)));
    end
    if EM(1,1) > 0
        disp('L_s2: '+string(EM(3,1))+' | '+string(EM(1,1)));
    elseif EM(1,1) < 0
        disp('C_s2: '+string(EM(3,1))+' | '+string(EM(1,1)));
    end
end
%--------------------------------------------------function-------------------------------------------------
%matching network sweep calcualtion
function Zout = sweep_matching(Im,Re,f,EM)
    R_coil = 0;
    w = 2*pi*f;
    %Load
    ZL = complex(Re,Im);
    %zero series element
    s0 = 0;
    if EM(1,5) == 0
        s0 = complex(0,0);
    elseif EM(1,5) > 0
        s0 = complex(R_coil,w*EM(3,5));
    elseif EM(1,5) < 0
        s0 = complex(0,-1/(w*EM(3,5)));
    end
    %calculate first parallel element
    p1 = 0;
    if EM(2,4) == 0
        p1 = 0;
    elseif EM(2,4) > 0
        p1 = complex(R_coil,w*EM(4,4));
    elseif EM(2,4) < 0
        p1 = complex(0,-1/(w*EM(4,4)));
    end
    %calculate first series element
    s1 = 0;
    if EM(1,3) == 0
        s1 = complex(0,0);
    elseif EM(1,3) > 0
        s1 = complex(R_coil,w*EM(3,3));
    elseif EM(1,3) < 0
        s1 = complex(0,-1/(w*EM(3,3)));
    end
    %calcualte second parallel element
    p2 = 0;
    if EM(2,2) == 0
        p2 = 0;
    elseif EM(2,2) > 0
        p2 = complex(R_coil,w*EM(4,2));
    elseif EM(2,2) < 0
        p2 = complex(0,-1/(w*EM(4,2)));
    end
    %calucalte second series element
    s2 = 0;
    if EM(1,1) == 0
        s2 = complex(0,0);
    elseif EM(1,1) > 0
        s2 = complex(R_coil,w*EM(3,1));
    elseif EM(1,1) < 0
        s2 = complex(0,-1/(w*EM(3,1)));
    end
    %calucalte the the matched frequency dependent impedance
    Z1 = ZL+s0;
    if p1 == 0
        Z2 = Z1+s1;
    else
        Z2 = (p1*Z1)/(p1+Z1)+s1;
    end
    if p2 == 0
        Zout = Z2+s2;
    else
        Zout = (p2*Z2)/(p2+Z2)+s2;
    end
end
%--------------------------------------------------function-------------------------------------------------
%matching network component calculation
function Zout = match_Z(Im,Re,f0,filter)
    EM = zeros(6,6);
    w = 2*pi*f0;
    EM(1,6) = Im;
    %EM(1,6) = -15;
    EM(2,6) = Re;
    %EM(2,6) = 10;
    EM(3,6) = -1/(w*EM(1,6));
    EM(4,6) = EM(2,6);
    %convert series to parallel
    Q0p = 1/(w*EM(3,6)*EM(4,6));
    R0p = EM(4,6)*(1+Q0p^2);
    C0p = EM(3,6)/(1+1/Q0p^2);
    %find suitable filter
    if filter=="resS"
        EM(1,5) = -EM(1,6);
        EM(3,5) = EM(1,5)/w;
        EM(5,1) = complex(EM(2,6),EM(1,6))+complex(0,EM(1,5));
    elseif filter=="resP"
        EM(2,4) = 1/(w*C0p);
        EM(4,4) = EM(2,4)/w;
        sL = complex(EM(2,6),EM(1,6));
        p1 = complex(0,EM(2,4));
        EM(5,1) = (sL*p1)/(sL+p1);
    elseif filter=="LS"
        if EM(4,6) > 50
            EM(1,5) = -EM(1,6);
            EM(3,5) = EM(1,5)/w;
            EM(6,3) = (EM(4,6)/50-1)^0.5;
            EM(2,4) = EM(4,6)/EM(6,3);
            EM(4,4) = EM(2,4)/w;
            EM(1,3) = -EM(6,3)*50;
            EM(3,3) = -1/(w*EM(1,3));
            sL = complex(EM(2,6),EM(1,6));
            s0 = complex(0,EM(1,5));
            p1 = complex(0,EM(2,4));
            s1 = complex(0,EM(1,3));
            EM(5,1) = ((sL+s0)*p1)/((sL+s0)+p1)+s1;
        elseif EM(4,6) <= 50
            EM(6,3) = (50/EM(4,6)-1)^0.5;
            s0 = -EM(6,3)*EM(4,6);
            EM(1,5) = s0-EM(1,6);
            if EM(1,5) > 0
                EM(3,5) = EM(1,5)/w;
            else
                EM(3,5) = -1/(w*EM(1,5));
            end
            EM(2,4) = 50/EM(6,3);
            EM(4,4) = EM(2,4)/w;
            sL = complex(EM(2,6),EM(1,6));
            s1 = complex(0,EM(1,5));
            p1 = complex(0,EM(2,4));
            EM(5,1) = ((sL+s1)*p1)/((sL+s1)+p1);
        end
    elseif filter =="LP"
        if R0p > 50
            EM(6,3) = (R0p/50-1)^0.5;
            So = -R0p/EM(6,3);
            Lo = -1/(w*C0p);
            EM(2,4) = (Lo*So)/(Lo-So);
            if EM(2,4) > 0
                EM(4,4) = EM(2,4)/w;
            else
                EM(4,4) = -1/(w*EM(2,4));
            end
            EM(1,3) = EM(6,3)*50;
            EM(3,3) = EM(1,3)/w;
            sL = complex(EM(2,6),EM(1,6));
            p1 = complex(0,EM(2,4));
            s1 = complex(0,EM(1,3));
            EM(5,1) = (sL*p1)/(sL+p1)+s1;
        elseif R0p <= 50
            EM(6,3) = (50/R0p-1)^0.5;
            EM(2,4) = 1/(w*C0p);
            EM(4,4) = EM(2,4)/w;
            EM(1,3) = EM(6,3)*R0p;
            EM(3,3) = EM(1,3)/w;
            EM(2,2) = -50/EM(6,3);
            EM(4,2) = -1/(w*EM(2,2));
            sL = complex(EM(2,6),EM(1,6));
            p0 = complex(0,EM(2,4));
            s1 = complex(0,EM(1,3));
            p1 = complex(0,EM(2,2));
            EM(5,4) = (sL*p0)/(sL+p0);
            EM(5,1) = ((EM(5,4)+s1)*p1)/((EM(5,4)+s1)+p1);
        end
    elseif filter=="LL"
        if EM(4,6) > 50
            Rv = (50*EM(4,6))^0.5;
            EM(6,3) = (EM(4,6)/Rv-1)^0.5;
            EM(6,1) = (Rv/50-1)^0.5;
            EM(1,5) = -EM(1,6);
            EM(3,5) = EM(1,5)/w;
            EM(2,4) = EM(4,6)/EM(6,3);
            EM(4,4) = EM(2,4)/w;
            EM(1,3) = -EM(6,3)*Rv;
            EM(3,3) = -1/(w*EM(1,3));
            EM(2,2) = -Rv/EM(6,1);
            EM(4,2) = -1/(w*EM(2,2));
            EM(1,1) = 50*EM(6,1);
            EM(3,1) = EM(1,1)/w;
            sL = complex(EM(2,6),EM(1,6));
            s0 = complex(0,EM(1,5));
            p1 = complex(0,EM(2,4));
            s1 = complex(0,EM(1,3));
            p2 = complex(0,EM(2,2));
            s2 = complex(0,EM(1,1));
            Z1 = ((sL+s0)*p1)/((sL+s0)+p1)+s1;
            EM(5,1) = (Z1*p2)/(Z1+p2)+s2;
        elseif EM(4,6) < 50
            Rv = (50*EM(4,6))^0.5;
            EM(6,3) = (Rv/EM(4,6)-1)^0.5;
            EM(6,1) = (50/Rv-1)^0.5;
            s0 = -EM(6,3)*EM(4,6);
            EM(1,5) = s0-EM(1,6);
            if EM(1,5) > 0
                EM(3,5) = EM(1,5)/w;
            else
                EM(3,5) = -1/(w*EM(1,5));
            end
            EM(2,4) = Rv/EM(6,3);
            EM(4,4) = EM(2,4)/w;
            EM(1,3) = EM(6,1)*Rv;
            EM(3,3) = EM(1,3)/w;
            EM(2,2) = -50/EM(6,1);
            EM(4,2) = -1/(w*EM(2,2));
            sL = complex(EM(2,6),EM(1,6));
            s1 = complex(0,EM(1,5));
            p1 = complex(0,EM(2,4));
            s2 = complex(0,EM(1,3));
            p2 = complex(0,EM(2,2));
            Z1 = ((sL+s1)*p1)/((sL+s1)+p1)+s2;
            EM(5,1) = (p2*Z1)/(p2+Z1);
        end
    elseif filter=="LLS"
        EM(6,5) = -EM(1,6)/EM(4,6);
        Rv0 = EM(4,6)*(1+EM(6,5)^2);
        if Rv0 > 50
            Rv = (50*Rv0)^0.5;
            EM(6,3) = (Rv0/Rv-1)^0.5;
            EM(6,1) = (Rv/50-1)^0.5;
            BL = EM(6,5)/Rv0;
            BC = -EM(6,3)/Rv0;
            EM(2,4) = 1/(BL+BC);
            if EM(2,4) > 0
                EM(4,4) = EM(2,4)/w;
            else
                EM(4,4) = -1/(w*EM(2,4));
            end
            EM(1,3) = EM(6,3)*Rv;
            EM(3,3) = EM(1,3)/w;
            EM(2,2) = Rv/EM(6,1);
            EM(4,2) = EM(2,2)/w;
            EM(1,1) = -EM(6,1)*50;
            EM(3,1) = -1/(w*EM(1,1));
            sL = complex(EM(2,6),EM(1,6));
            p1 = complex(0,EM(2,4));
            s1 = complex(0,EM(1,3));
            p2 = complex(0,EM(2,2));
            s2 = complex(0,EM(1,1));
            Z1 = (sL*p1)/(sL+p1)+s1;
            EM(5,1) = (Z1*p2)/(Z1+p2)+s2;
        elseif Rv0 <50
            EM(6,3) = -EM(1,6)/EM(4,6);
            Rv = EM(4,6)*(1+EM(6,3)^2);
            if Rv > 50
                EM(6,1) = (Rv/50-1)^0.5;
                BL = EM(6,3)/Rv;
                BC = -EM(6,1)/Rv;
                EM(2,4) = 1/(BL+BC);
                if EM(2,4) > 0
                    EM(4,4) = EM(2,4)/w;
                else
                    EM(4,4) = -1/(w*EM(2,4));
                end
                EM(1,3) = EM(6,1)*50;
                EM(3,3) = EM(1,3)/w;
                sL = complex(EM(2,6),EM(1,6));
                p1 = complex(0,EM(2,4));
                s1 = complex(0,EM(1,3));
                EM(5,1) = (sL*p1)/(sL+p1)+s1;
            elseif Rv < 50
                EM(6,1) = (50/Rv-1)^0.5;
                EM(2,4) = Rv/EM(6,3);
                EM(4,4) = EM(2,4)/w;
                EM(1,3) = EM(6,1)*Rv;
                EM(3,3) = EM(1,3)/w;
                EM(2,2) = -50/EM(6,1);
                EM(4,2) = -1/(w*EM(2,2));
                sL = complex(EM(2,6),EM(1,6));
                p1 = complex(0,EM(2,4));
                s2 = complex(0,EM(1,3));
                p2 = complex(0,EM(2,2));
                Z1 = (sL*p1)/(sL+p1)+s2;
                EM(5,1) = (Z1*p2)/(Z1+p2);
            end
        end
    end
    Zout = EM;
end
%--------------------------------------------------function-------------------------------------------------
%smoothen and convert calculated data
function Zout = smooth_theo_Z(Zin,amount)
    Z_imag = smooth(smooth(imag(Zin),amount,'rloess'),4);
    Z_real = (real(Zin).*4)';
    Zout = complex(Z_real,Z_imag);
end
%--------------------------------------------------function-------------------------------------------------
%calculate the Impedance and reflection from outside towards center
function Zin = calculate_Z(freq,Zi,US,US_shape,name,piezo_pos,C0,N)
    for f = 1:piezo_pos
        if f<=piezo_pos
            if f==1
                Zi(2,f) = Zi(1,f);
            end
            [Zi(2,f+1),Zi(3,f)] = trans(Zi(1,f+1),Zi(2,f),wave_const(freq,US.(string(name(f))).v),US.(string(name(f))).t);
            %%%%%%%%only for visulation%%%%%%%%%
            [a,Zi(3,f)] = trans(Zi(1,f+1),Zi(1,f),wave_const(freq,US.(string(name(f))).v),US.(string(name(f))).t);
        end
    end
    %calculate right of piezo towards center
    for b = (US_shape(1,1)+4):-1:(piezo_pos+4)
        if b==US_shape(1,1)+4
            Zi(2,b) = Zi(1,b);
        end
        [Zi(2,b-1),Zi(3,b)] = trans(Zi(1,b-1),Zi(2,b),wave_const(freq,US.(string(name(b-4))).v),US.(string(name(b-4))).t);
        %%%%%%%%only for visulation%%%%%%%%%
        [a,Zi(3,b)] = trans(Zi(1,b-1),Zi(1,b),wave_const(freq,US.(string(name(b-4))).v),US.(string(name(b-4))).t);
    end
    %calulcate the center accoustic
    piezo_cent = piezo_pos+2;
    [ZT, ZS] = half_trans(Zi(1,piezo_cent-1),wave_const(freq,US.(string(name(piezo_pos))).v),US.(string(name(piezo_pos))).t);
    %combine the piezo impedance with the surrounding
    Zi(2,piezo_cent-1) = ZT+Zi(2,piezo_cent-2);
    Zi(2,piezo_cent+1) = ZT+Zi(2,piezo_cent+2);
    %Impedance on secondary side
    Zi(2,piezo_cent) = ZS+(Zi(2,piezo_cent-1)*Zi(2,piezo_cent+1))/(Zi(2,piezo_cent-1)+Zi(2,piezo_cent+1));
    %transform impedance to the primary side
    Ztrans = Zi(2,piezo_cent)*(1/N)^2;
    ZC = 1/(1i*2*pi*freq*C0);
    %input impedance
    Zi(1,piezo_cent) = (ZC*(-ZC+Ztrans))/(ZC+(-ZC+Ztrans));
    
    Zin=Zi;
end
%--------------------------------------------------function-------------------------------------------------
%calculate the center impedances of the piezo
function [ZT, ZS] = half_trans(Z0,beta,t)
    ZT = 1i*Z0*tan(beta*t/2);
    ZS = -1i*Z0*1/sin(beta*t);
end
%--------------------------------------------------function-------------------------------------------------
%calculate the impedance and reflection that are attached to the piezo
function [Zin,gamma] = trans(Z0,ZL,beta,t)
    Zin = Z0*(ZL+1i*Z0*tan(beta*t))/(Z0+1i*ZL*tan(beta*t));
    gamma = (ZL-Z0)/(ZL+Z0);
end
%--------------------------------------------------function-------------------------------------------------
%calculate the wave properagation
function beta = wave_const(f,v)
    beta = (2*pi*f)/v;
end