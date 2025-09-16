clear
ecoli = load('D:\OneDrive - smail.iitm.ac.in\SprintCore\Bigg models\e_coli_core.mat');
ecoli.model.c(25)=1;
ecoli.model.rev = ecoli.model.lb<0;

recon = load('D:\OneDrive - smail.iitm.ac.in\SprintCore\Bigg models\Recon3D.mat');
recon.model.c(2502)=1;
recon.model.rev = recon.model.lb<0;

ecoli_sol = optimizeCbModel(ecoli.model);
recon_sol = optimizeCbModel(recon.model);

ecoli_sol_1 = optimizeCbModel(ecoli.model,'max','one');
recon_sol_1 = optimizeCbModel(recon.model,'max','one');

% displaying number of non-zero fluxes in regular FBA
disp(['Number of non-zero fluxes in E. coli model: ', num2str(nnz(ecoli_sol.x))]);
disp(['Number of non-zero fluxes in Recon3D model: ', num2str(nnz(recon_sol.x))]);

% displaying the objective value
disp(['Objective value for E. coli model: ', num2str(ecoli_sol.f)]);
disp(['Objective value for Recon3D model: ', num2str(recon_sol.f)]);

% displaying the number of non-zero fluxes in minNorm FBA
disp(['Number of non-zero fluxes in E. coli model (minNorm): ', num2str(nnz(ecoli_sol_1.x))]);
disp(['Number of non-zero fluxes in Recon3D model (minNorm): ', num2str(nnz(recon_sol_1.x))]);

recon.model.lb(2502) = recon_sol.f;
ecoli.model.lb(25) = ecoli_sol.f;

% running the flux reducer
recon_flux_reduced = FluxReducer(recon.model, 1);
ecoli_flux_reduced = FluxReducer(ecoli.model, 1);

% displaying the number of non-zero fluxes after reduction
disp(['Number of non-zero fluxes in E. coli model after reduction: ', num2str(nnz(ecoli_flux_reduced))]);
disp(['Number of non-zero fluxes in Recon3D model after reduction: ', num2str(nnz(recon_flux_reduced))]);

% displaying the reduced objective values
disp(['Reduced objective value for E. coli model: ', num2str(ecoli_flux_reduced(25))]);
disp(['Reduced objective value for Recon3D model: ', num2str(recon_flux_reduced(2502))]);