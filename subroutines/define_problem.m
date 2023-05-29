function problem = define_problem(varargin)

if nargin == 0
    name = 'doubleglazing';
end

switch varargin{1}
    case 'doubleglazing'
        problem.linear = 1;
        problem.wind_fn = 'recirculating-uncertain';
        problem.nWind = 4;
        problem.sigmaWind = 0.5;

        problem.diff_fn = 'fixed';
        problem.viscosity = 1e-1;        
        problem.nDiff = 0;

        problem.bc = 'hotwall';
        problem.tau = 1e-1;
        problem.f = 'zero';
    case 'doubleglazing-vardiff'
        problem.linear = 1;
        problem.wind_fn = 'recirculating-uncertain';
        problem.nWind = 4;
        problem.sigmaWind = 0.5;

        problem.diff_fn = 'scalar-var';
        problem.viscosity = 1e-1;
        problem.viscosity_var = 9e-2;
        problem.nDiff = 1;

        problem.bc = 'hotwall';
        problem.tau = 1e-1;
        problem.f = 'zero';
    case 'doubleglazing-timebc'
        problem.linear = 1;
        problem.wind_fn = 'recirculating-uncertain';
        problem.nWind = 4;
        problem.sigmaWind = 0.5;

        problem.diff_fn = 'scalar-var';
        problem.viscosity = 1e-1;
        problem.viscosity_var = 9e-2;
        problem.nDiff = 1;

        problem.bc = 'hotwall-timebc';
        problem.tau = 1e-1;
        problem.f = 'zero';
    case 'doubleglazing-16'
        problem.linear = 1;
        problem.wind_fn = 'recirculating-uncertain-decay';
        problem.nWind = 16;
        problem.sigmaWind = 0.5;

        problem.diff_fn = 'fixed';
        problem.viscosity = 1e-1;
        problem.nDiff = 0;

        problem.bc = 'hotwall';
        problem.tau = 1e-1;
        problem.f = 'zero';
    case 'doubleglazing-20'
        problem.linear = 1;
        problem.wind_fn = 'jomp';
        problem.nWind = 0.5;
        problem.sigmaWind = 0.5;

        problem.diff_fn = 'fixed';
        problem.viscosity = 1e-1;
        problem.nDiff = 0;

        problem.bc = 'hotwall';
        problem.tau = 1e-1;
        problem.f = 'zero';
    case 'doubleglazing-32'
        problem.linear = 1;
        problem.wind_fn = 'jomp32';
        problem.nWind = 32;
        problem.sigmaWind = 1;

        problem.diff_fn = 'fixed';
        problem.viscosity = 1e-1;
        problem.nDiff = 0;

        problem.bc = 'hotwall';
        problem.tau = 1e-1;
        problem.f = 'zero';
    case 'doubleglazing-64'
        problem.linear = 1;
        problem.wind_fn = 'jomp64';
        problem.nWind = 64;
        problem.sigmaWind = 1;

        problem.diff_fn = 'fixed';
        problem.viscosity = 1e-1;
        problem.nDiff = 0;

        problem.bc = 'hotwall';
        problem.tau = 1e-1;
        problem.f = 'zero';
    case 'doubleglazing-64-rf'
        problem.linear = 1;
        problem.wind_fn = 'djs-rf';
        problem.nWind = 64;
        problem.sigmaWind = 1;

        problem.diff_fn = 'fixed';
        problem.viscosity = 1e-1;
        problem.nDiff = 0;

        problem.bc = 'hotwall';
        problem.tau = 1e-1;
        problem.f = 'zero';
    case 'doubleglazing-no-rv'
        problem.linear = 1;
        problem.wind_fn = 'recirculating';
        problem.nWind = 0;
        problem.sigmaWind = 0;

        problem.diff_fn = 'fixed';
        problem.viscosity = 1e-1;
        problem.nDiff = 1;

        problem.bc = 'hotwall';
        problem.tau = 1e-1;
        problem.f = 'zero';

    case 'vidlickova'
        problem.linear = 1;
        problem.wind_fn = 'zero';
        problem.nWind = 0;

        problem.diff_fn = 'vidlickova';
        problem.nDiff = 3;

        problem.bc = 'dirichlet';
        problem.f = 'vidlickova';

    case 'eigel'
        problem.linear = 1;
        problem.wind_fn = 'zero';
        problem.nWind = 0;

        problem.diff_fn = 'eigel';
        problem.nDiff = 10;

        problem.bc = 'dirichlet';
        problem.f = 'one';

    case 'xwind'
        problem.linear = 1;
        problem.wind_fn = 'xonly';
        problem.nWind = 0;
        problem.windmag = 2;

        problem.diff_fn = 'eigel';
        problem.nDiff = 10;

        problem.bc = 'dirichlet';
        problem.f = 'square';
        problem.fwidth = 0.1;
        problem.fmag = 20;

    case 'groundwater'
        problem.linear = 1;
        problem.wind_fn = 'zero'; %'recirculating' ; %'xonly-weak'; % 'recirculating-uncertain'; %'pouiseuille'; %% 'prevailing-perturbed'; %'recirculating-uncertain'; %'zero'; %'recirculating'; %'xonly-weak'; %  'zero'; %%-perturbed'; %'recirculating-uncertain'; % 'recirculating'; %'recirculating-alt'; %'recirculating-uncertain'; %'zero';'zero'; %
        problem.nWind = 0;
        problem.sigmaWind = 0.5;

        problem.diff_fn = 'eigel'; %'fixed'; %'complex'; %%'fixed'; %'vidlickova'; %'eigel';
        problem.viscosity = 1e-1;

        problem.nDiff = 10;

        problem.linear =0;
        problem.wind_fn = 'groundwater';
        problem.diff_fn = 'groundwater';
        problem.n = 5;

        problem.bc = 'dirichlet'; %'pipe'; %'neumann'; %'hotwall'; %%
        problem.tau = 1e-1;

        problem.f = 'one'; %'square'; % 'zero'; %'vidlickova'; %'one'; %'square'; %'circle'
    otherwise
        error('No such problem exits');
end

problem.n = problem.nWind+problem.nDiff;

