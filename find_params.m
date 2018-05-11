store_params = load('store_params.mat');
store_acc_val = load('store_acc_val.mat');
store_params = struct2cell(store_params);
store_params = cell2mat(store_params);
store_acc_val = struct2cell(store_acc_val);
store_acc_val = cell2mat(store_acc_val);

% store_params is  2000 row x 7 column x 8 back matrix of all tried
% parameters (ksr1, krs2, gs1, gr1), the corresponding Nsens and Nres they 
% out put, and the error between that and the experimental value
 
% store_acc_val is a  8 row by 7 column matrix where each row is
%ksr1, krs2, gs1, gr1, Nsens_model, Nres_model, error in guess

Nsens_model =vertcat(Nsens(1), store_acc_val(:,5));
Nres_model = vertcat(Nres(1),store_acc_val(:,6));

ksr_model = store_acc_val(:,1);
krs_model = store_acc_val(:,2);
gs_model = store_acc_val(:,3);
gr_model = store_acc_val(:,4);

eff_ksr = ksr_model-krs_model;

figure(2)
set(gca,'LineWidth',1.5,'FontSize',12);
hold off
plot(1:1:8, ksr_model,'-g', 'LineWidth',2)
hold on
plot(1:1:8, krs_model,'-r', 'LineWidth',2)
plot(1:1:8, eff_ksr, '-k','LineWidth',2)
plot(1:1:8, gs_model,'-y', 'LineWidth',2)
plot(1:1:8, gr_model, '-b','LineWidth',2)
xlabel('Time Post Treatment (Weeks)', 'FontSize', 14)
ylabel('Parameter Value', 'FontSize', 14)
title('Final Accepted Parameter Values at Each Interval', 'FontSize', 14)
%legend('transition s->r', 'transition r->s', 's->r effective', 'sensitive net growth', 'resistant net growth')

figure (1)
set(gca,'LineWidth',1.5,'FontSize',18);
hold off
plot(0:1:8, log(Nsens_model),'-y', 'LineWidth',1)
hold on
plot(0:1:8, log(Nres_model), '-b', 'LineWidth',1)
plot(0:1:8, log(Nsens), 'y*', 'LineWidth',2)
plot(0:1:8, log(Nres), 'b*', 'LineWidth',2)
xlabel('Time Post Treatment(Weeks)', 'FontSize', 14)
ylabel( 'Log Tumor Cell Number', 'FontSize',14)
%legend('Sensitive Model', 'Resistant Model', 'Sensitive Expt', 'Resistant Expt')
title('Model Fit versus Experimentall Derived Estimate', 'FontSize',14)
%%
dt = .1;
alpha = 100;
for i = 1:8
    for n = 1:2000
        if i == 1
        Err(n,1,i) = abs(store_params(n,7, i)) +alpha*abs( (1/dt*(store_params(n,1, i+1)-store_params(n,1,i))+...
        1/dt*(store_params(n,2, i+1)-store_params(n,2,i))));
%         1/dt*(store_params(n,3, i+1)-store_params(n,3,i))+...
%         1/dt*(store_params(n,4, i+1)-store_params(n,4,i))));
        elseif i == 8
        Err(n,1,i) = store_params(n,7, i) + alpha*abs((1/dt*((store_params(n,1, i)-store_params(n,1,i-1)))+...
        1/dt*(store_params(n,2, i)-store_params(n,2,i-1))));
%         1/dt*(store_params(n,3, i)-store_params(n,3,i-1))+...
%         1/dt*(store_params(n,4, i)-store_params(n,4,i-1))));
        else
    Err(n,1,i) = store_params(n,7, i) + alpha*abs((1/dt*((store_params(n,1, i+1)-store_params(n,1,i-1)))+...
        1/dt*(store_params(n,2, i+1)-store_params(n,2,i-1))));
%         1/dt*(store_params(n,3, i+1)-store_params(n,3,i-1))+...
%         1/dt*(store_params(n,4, i+1)-store_params(n,4,i-1))));
        end

    end

end
%%
% Add new error to store params, then find minimum
for i = 1:8
    min_err(i,1) = min((Err(:,1,i)));
    % need to find index that this occurs ** I think this is wrong
    % currently
    min_err(i,2) = find(Err == min(Err(:,1, i)));
    min_err(i,3) = min_err(i,2) - (i-1)*2000;
end

% Now you have indices, and want to extract them and make into a new
% vector to store

for i = 1:8
    store_new_params(i,:) = store_params(min_err(i,3), :, i);
end

ksr_model2 = store_new_params(:,1);
krs_model2 = store_new_params(:,2);
gs_model2 = store_new_params(:,3);
gr_model2 = store_new_params(:,4);

eff_ksr2 = ksr_model2-krs_model2;

% Need to find a better way 

figure(3)
set(gca,'LineWidth',1.5,'FontSize',12);
hold off
plot(1:1:8, ksr_model2,'-g', 'LineWidth',2)
hold on
plot(1:1:8, krs_model2,'-r', 'LineWidth',2)
plot(1:1:8, eff_ksr2, '-k','LineWidth',2)
plot(1:1:8, gs_model2,'-y', 'LineWidth',2)
plot(1:1:8, gr_model2, '-b','LineWidth',2)
xlabel('Time Post Treatment (Weeks)')
ylabel('Parameter Value')
title('Final Accepted Parameter Values at Each Interval')
legend('transition s->r', 'transition r->s', 's->r effective', 'sensitive net growth', 'resistant net growth')

Tot_err = sum(store_acc_val(:,7))