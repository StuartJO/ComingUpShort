%% Description:
% This script performs random smoothing of values across Schaefer 400 cortical parcellations 
% on both hemispheres of the brain using fsaverage surface topology data.
% It:
% 1. Loads left hemisphere (LH) surface data and parcellations.
% 2. Computes neighboring parcels for each LH parcel.
% 3. Performs iterative smoothing of random values across parcels based on neighbors.
% 4. Projects these smoothed LH values onto the right hemisphere (RH) parcellation
% 5. Calculates the correlation matrix between smoothed LH and RH random values.
% 

%% Load fsaverage surface data (faces, vertices, and Schaefer 400 parcellation)
load('fsaverage_surface_data.mat');   % Contains lh_faces, lh_verts, and lh_scha400 parcellation data

faces = lh_faces;                      % Face connectivity for LH surface
verts = lh_verts;                      % Vertex coordinates for LH surface
parc_id = Scha7_parcs.lh_scha400;      % Parcel ID for each vertex (0 if not assigned)

nParcs = max(parc_id);                 % Number of parcels (LH)

iters = 5;                             % Number of smoothing iterations
nFeatures = 20;                        % Number of features to create
ParcNei = cell(1, nParcs);             % Preallocate cell array for parcel neighbors

RAND_SMOOTH_L = zeros(nFeatures, 200);        % Preallocate matrix to store LH smoothed values (20 random instances, 200 parcels)

%% Smooth random values across LH parcels
for F = 1:nFeatures
    current_smoothed = rand(1, 200);   % Generate random values for each LH parcel
    new_smoothed = current_smoothed;

    for iter = 1:iters
        % Compute smoothing by averaging parcel values with their neighbors
        for i = 1:nParcs
            if isempty(ParcNei{i})
                % Find neighbors if not already computed
                parc_verts = find(parc_id == i);
                I = sum(ismember(faces, parc_verts), 2) > 0;  
                Uverts = unique(faces(I, :));                 
                Uparcs = unique(parc_id(Uverts));             
                Uparcs(Uparcs == i | Uparcs == 0) = [];       
                ParcNei{i} = Uparcs;                          
            end
            new_smoothed(i) = mean(current_smoothed(ParcNei{i})); 
        end
        current_smoothed = new_smoothed; % Update for next iteration
    end

    RAND_SMOOTH_L(F, :) = current_smoothed; % Store result
end

%% Load fsLR surface files for both hemispheres
% fsLR has vertex correspondence across hemispheres, while fsaverage does
% not. Therefore we can use fsLR to project one hemipshere onto the other
G = gifti('S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii');
lh_verts = double(G.vertices);
lh_faces = double(G.faces);

G = gifti('S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii');
rh_verts = double(G.vertices);
rh_faces = double(G.faces);

%% Load fsLR Schaefer 400 parcellation for both hemispheres
lh_scha400_fsLR = dlmread('fsLR_32k_Schaefer400-lh.txt');       
rh_scha400_fsLR = dlmread('fsLR_32k_Schaefer400-rh.txt') - 200; 
rh_scha400_fsLR(rh_scha400_fsLR == -200) = 0;                   % Set unassigned vertices to zero

RAND_SMOOTH_R = zeros(nFeatures, 200); % Preallocate smoothed RH matrix

%% Map LH smoothed values to RH parcels and compute average value for each RH parcel
for F = 1:nFeatures
    % Assign LH smoothed values to LH vertices and project onto RH
    rh_scha400_fsLR_random_vals = my_changem(lh_scha400_fsLR, RAND_SMOOTH_L(F, :), 1:200);

    % Average projected values within each RH parcel
    for i = 1:200
        RAND_SMOOTH_R(F, i) = mean(rh_scha400_fsLR_random_vals(rh_scha400_fsLR == i));
    end
end

%% Compute correlation between LH and RH smoothed random values
C = corr([RAND_SMOOTH_L, RAND_SMOOTH_R], 'Type', 'Pearson');
