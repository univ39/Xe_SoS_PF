function phase = generateSmoothPhase(sz_or_mask, varargin)
% generateSmoothPhase  Create a smoothly varying 3D phase map (radians)
%
% Usage:
%   phase = generateSmoothPhase([nx,ny,nz])
%   phase = generateSmoothPhase(nx,ny,nz)
%   phase = generateSmoothPhase(mask)        % mask is logical or numeric 3D
%
% Name-value options:
%   'PhaseAmp'      - max absolute phase (radians). Default pi/2.
%   'RandPhaseAmp'  - amplitude of low-frequency random component. Default 0.4.
%   'BumpAmp'       - gaussian bump amplitude. Default 0.7.
%   'BumpSigma'     - bump width in normalized coords (0..1). Default 0.25.
%   'PlaneCoef'     - 1x3 vector [ax ay az] plane coefficients. Default [0.6 -0.4 0.2].
%   'QuadScale'     - scale for quadratic term. Default 0.9.
%   'Seed'          - random seed (numeric) or [] for shuffle. Default [].
%   'ApplyMask'     - if true and first arg was a mask, phase outside mask set to 0. Default true.
%
% Output:
%   phase  - nx-by-ny-by-nz real array (radians)
%
% Example:
%   % create phase for 256x256x64
%   ph = generateSmoothPhase([256,256,64],'PhaseAmp',pi,'Seed',123);
%   imagesc(angle(exp(1i*ph(:,:,32)))); colorbar; title('phase slice');
%
%   % create phase using an existing mask (apply inside mask only)
%   mask = (rand(128,128,32) > 0.6);
%   ph = generateSmoothPhase(mask,'ApplyMask',true);
%
% Author: ChatGPT
% Date:   2025-10-24

%% parse inputs
p = inputParser;
p.KeepUnmatched = true;

validScalar = @(x) isnumeric(x) && isscalar(x);
addParameter(p,'PhaseAmp',pi/2,@isnumeric);
addParameter(p,'RandPhaseAmp',0.4,@isnumeric);
addParameter(p,'BumpAmp',0.7,@isnumeric);
addParameter(p,'BumpSigma',0.25,@isnumeric);
addParameter(p,'PlaneCoef',[0.6 -0.4 0.2],@(x)isnumeric(x)&&numel(x)==3);
addParameter(p,'QuadScale',0.9,@isnumeric);
addParameter(p,'Seed',[],@(x)(isempty(x) || (isnumeric(x)&&isscalar(x))));
addParameter(p,'ApplyMask',true,@islogical);

% handle flexible first arg (size vector, three scalars, or mask)
if isvector(sz_or_mask) && numel(sz_or_mask)==3 && all(sz_or_mask==round(sz_or_mask))
    nx = sz_or_mask(1); ny = sz_or_mask(2); nz = sz_or_mask(3);
    maskProvided = false;
    providedMask = [];
elseif isnumeric(sz_or_mask) && isscalar(sz_or_mask) && nargin>=3
    nx = sz_or_mask; ny = varargin{1}; nz = varargin{2};
    % remove the two scalars from varargin for parser
    varargin(1:2) = [];
    maskProvided = false;
    providedMask = [];
elseif (islogical(sz_or_mask) || isnumeric(sz_or_mask)) && ndims(sz_or_mask)==3
    providedMask = sz_or_mask ~= 0;
    [nx, ny, nz] = size(providedMask);
    maskProvided = true;
else
    error('First argument must be [nx,ny,nz], nx,ny,nz, or a 3D mask/volume.');
end

parse(p, varargin{:});
opts = p.Results;

% RNG
if ~isempty(opts.Seed)
    rng(opts.Seed);
else
    rng('shuffle');
end

%% coordinate grids in normalized [-1,1]
[xg, yg, zg] = ndgrid(linspace(-1,1,nx), linspace(-1,1,ny), linspace(-1,1,nz));

%% deterministic components
% plane
plane = opts.PlaneCoef(1)*xg + opts.PlaneCoef(2)*yg + opts.PlaneCoef(3)*zg;

% quadratic term (spherical-ish)
quad = opts.QuadScale * (xg.^2 + yg.^2 + 0.5*zg.^2);

% gaussian bump (off-center for visual interest)
xc = 0.3; yc = -0.2; zc = 0.0;
bump = opts.BumpAmp * exp( - ( (xg - xc).^2 + (yg - yc).^2 + (zg - zc).^2 ) / (2*opts.BumpSigma^2) );

%% low-frequency random field (smoothed Gaussian noise)
R = randn(nx,ny,nz);

% smoothing widths (fractions of dims -> pixels)
sigma_frac = 0.12;            % fraction of each dimension for smoothing
sigma_px = max(1, round(sigma_frac .* [nx, ny, nz]));

gx = gauss1d(sigma_px(1));
gy = gauss1d(sigma_px(2));
gz = gauss1d(sigma_px(3));

% separable smoothing
Rsm = convn(R, reshape(gx,[],1,1), 'same');
Rsm = convn(Rsm, reshape(gy,1,[],1), 'same');
Rsm = convn(Rsm, reshape(gz,1,1,[]), 'same');

% normalize to unit std then scale
if std(Rsm(:))>0
    Rsm = Rsm ./ std(Rsm(:));
end
randPhase = opts.RandPhaseAmp * Rsm;

%% combine
rawPhase = plane + quad + bump + randPhase;

% scale so maximum absolute value becomes PhaseAmp
if opts.PhaseAmp > 0
    maxabs = max(abs(rawPhase(:)));
    if maxabs > 0
        phase = rawPhase / maxabs * opts.PhaseAmp;
    else
        phase = rawPhase;
    end
else
    phase = rawPhase;
end

%% apply mask option (if user provided a mask)
if maskProvided && opts.ApplyMask
    phase(~providedMask) = 0;
end

end

%% helper: 1D gaussian kernel
function g = gauss1d(sigma)
    if sigma <= 0
        g = 1;
        return;
    end
    truncate = 4;
    halfw = ceil(truncate * sigma);
    x = -halfw:halfw;
    g = exp(-(x.^2)/(2*sigma^2));
    g = g / sum(g(:));
end
