function [D,new_T1_record,new_T2_record,new_B1_record] = generate_multimapping_dictionary(...
    T1_record,T2_record,B1_record,xDoc,userOutput)
% GENERATE_MULTIMAPPING_DICTIONARY: Generate a dictionary of multimapping
% signals. Each column (entry) of the dictionary is a multimapping signal
% with one of the given T1 values [ms], T2 values [ms] and B1 values, while
% each row of the dictionary corresponds to one of the image sampling times.

n_images = 10;  % Number of the multimapping images (image sampling times)
nr_startups = 5;  % Number of the excitations prior to the steady state

% Define the RF pulses (the actual imaging/readout flip angles will be defined later during B1 determination)
t2prep_rf = [90 180 90];
inv_rf = 180;
fa = xDoc.Root.Seq.Basic.FlipAngle.Value;
n_ex = xDoc.Root.Seq.KSpace.Segments.Value;

% Flip angles of the RF pulses prior to the steady state
SU_fa_list = zeros(1,nr_startups);
for i = 1:nr_startups
    SU_fa_list(i) = fa*i/nr_startups;
end

npulse = n_ex+nr_startups+4;  % Maximum number of the EPG entries for a cardiac cycle. 4 = nr RF in T2 prep + trigger delay
TR = xDoc.Root.Seq.Basic.TR.Value/1000;

% Define the phases for the RF pulses
phi_bSSFP = RF_phase_cycle(npulse,'balanced');
for i = 1:n_images
    phi(i,:) = phi_bSSFP;
end
phi_t2prep = d2r([0 90 180]);

% The T2 preps are performed first in cardiac cyle 8, 9 & 10
phi(8,1:3) = phi_t2prep;
phi(9,1:3) = phi_t2prep;
phi(10,1:3) = phi_t2prep;

% Determine the delays
% Delay_dur = RR interval minus all fixed durations
% TNT_dur = duration from last RF pulse in current cardiac cycle to the first
% one in the next cycle. In this case, no RF pulses are performed during
% inversion or T2 prep.
centersegment = 41;
if xDoc.Root.Seq.PPA.Method.Value == 1
    Index = find (userOutput(:,3) == userOutput(1,3));
    for i = 2:n_images
        TNT_dur(i-1) = userOutput(Index(i),6)-userOutput(Index(i)-1,6);
    end
    TNT_dur = TNT_dur/10-nr_startups*TR;
    TNT_dur(n_images) = 0;
else
    TNT_dur(1) = userOutput(1,5)+seq.TIdiff*10;
    for i = 1:n_images-1
        TNT_dur(i-1) = userOutput(seq.nview+1,6)-userOutput(seq.nview,6);
    end
    TNT_dur = TNT_dur/10;
end
TFEPP_dur = 20;
TI = cell2mat(xDoc.Root.Seq.Basic.TI.Value)/1000;
TFEPP_delay_dur(1) = TI(1)-(nr_startups+centersegment)*TR-TFEPP_dur*0.5;
TFEPP_delay_dur(5) = TI(2)-(nr_startups+centersegment)*TR-TFEPP_dur*0.5;

T2prep = cell2mat(xDoc.Root.Seq.PrepPulse.T2PrepDuration.Value);
TNT_dur(4) = TNT_dur(4)-TFEPP_dur*0.5-TFEPP_delay_dur(5);
TNT_dur(7) = TNT_dur(7)-T2prep(1);
TNT_dur(8) = TNT_dur(8)-T2prep(2);
TNT_dur(9) = TNT_dur(9)-T2prep(3);
TNT_dur(10) = 0;

PP_delay = zeros(n_images,2);
PP_delay(1,1) = TFEPP_dur*0.5+TFEPP_delay_dur(1);
PP_delay(5,1) = TFEPP_dur*0.5+TFEPP_delay_dur(5);
PP_delay(8,2) = T2prep(1);
PP_delay(9,2) = T2prep(2);
PP_delay(10,2) = T2prep(3);

% Calculate index of k0 acquisition of pulse train
k0_ind1 = nr_startups+centersegment;
k0_ind(1) = 1+k0_ind1;
k0_ind(2) = k0_ind1;
k0_ind(3) = k0_ind1;
k0_ind(4) = k0_ind1;
k0_ind(5) = 1+k0_ind1;
k0_ind(6) = k0_ind1;
k0_ind(7) = k0_ind1;
k0_ind(8) = 3+k0_ind1;
k0_ind(9) = 3+k0_ind1;
k0_ind(10) = 3+k0_ind1;

% Generate the dictionary
T1_record = T1_record(:);
T2_record = T2_record(:);
B1_record = B1_record(:);

T1_num = numel(T1_record);
T2_num = numel(T2_record);
B1_num = numel(B1_record);

para_list = zeros(T1_num*T2_num*B1_num,3);
para_num = 0;
for i = 1:T1_num
    for j = 1:T2_num
        for k = 1:B1_num

            if T1_record(i) > T2_record(j)
                para_num = para_num+1;
                para_list(para_num,1) = T1_record(i);
                para_list(para_num,2) = T2_record(j);
                para_list(para_num,3) = B1_record(k);
            end

        end
    end
end
para_list = para_list(1:para_num,:);

fprintf('Generate the dictionary...\n')

tic
D = zeros(n_images,para_num);
parfor para_ind = 1:para_num

    sim_im_val = zeros(n_images,1);

    % Scale the nominal flip angle (fa) with B1 factor in para_list(para_ind,3)
    fa2 = fa*para_list(para_ind,3);
    ss = fa2*ones([1 n_ex]);
    SU_fa_list2 = SU_fa_list*para_list(para_ind,3);

    ss = [SU_fa_list2 ss];

    alphas = [];

    alphas(8,:) = d2r([t2prep_rf ss 0]);
    alphas(9,:) = d2r([t2prep_rf ss 0]);
    alphas(10,:) = d2r([t2prep_rf ss 0]);

    alphas(1,:) = d2r([inv_rf ss 0 NaN NaN]);
    alphas(2,:) = d2r([ss 0 NaN NaN NaN]);
    alphas(3,:) = d2r([ss 0 NaN NaN NaN]);
    alphas(4,:) = d2r([ss 0 NaN NaN NaN]);
    alphas(5,:) = d2r([inv_rf ss 0 NaN NaN]);
    alphas(6,:) = d2r([ss 0 NaN NaN NaN]);
    alphas(7,:) = d2r([ss 0 NaN NaN NaN]);

    Mz = 1;

    % Loop over the cardiac cycles. Only keep track of Mz across the cardiac
    % cycles. Mxy is assumed to be spoiled between cycles
    for im_ind = 1:n_images

        n_rf = npulse;
        tmp = isnan(alphas(im_ind,:));
        if sum(tmp) > 0
            tmp2 = find(tmp == 1);
            n_rf = min(n_rf,min(tmp2)-1);
        end

        % Perform the EPG
        [~,Fn1,Zn1] = EPG_GRE_dict(alphas(im_ind,1:n_rf),phi(im_ind,1:n_rf),TR,para_list(para_ind,1),para_list(para_ind,2),TNT_dur(im_ind),PP_delay(im_ind,:),Mz,'kmax',inf);

        ssfp = ifftshift(size(Fn1,1)*(ifft(ifftshift(Fn1,1),[],1)),1);
        phiTR= linspace(-pi,pi,size(ssfp,1));
        [~,idx0] = min(abs(phiTR));
        Mxy = abs(ssfp(idx0,:));
        sign_Mz = sign(Zn1(1,:));
        Mxy = (Mxy.*sign_Mz);
        Mz = Zn1(1,end);

        sim_im_val(im_ind) = Mxy(k0_ind(im_ind));

    end
    D(:,para_ind) = sim_im_val;

end
toc

% Take the magnitude
D = abs(D);

% Normalization
D = D./sqrt(sum(D.^2,1));

new_T1_record = para_list(:,1);
new_T2_record = para_list(:,2);
new_B1_record = para_list(:,3);

end

