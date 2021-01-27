function JZiai_MSR10166_cd3(varargin) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JZIAI_MSR10166_CD3 processes a bright-field image from one slide. 
%
% Filename: JZiai_MSR10166_cd3.m
% Compile (if necessary, next two lines):
%   mcc -m JZiai_MSR10166_cd3 
% Use (within Matlab):  JZiai_MSR10166_cd3
% -inputdir 'NDP_Pcore2/Genentech_Webslide_server/MSRs/20170303_DDunlap_MSR10167/gslide_filename1.ndpi' 
%
% Output files:   
%   stats in .xls file
%
% Dependencies that must be in the same directory as this file.
%   -- none
%
% v1: new code.  from MSR228_cd3.  necroROI, subsets, xy coord for cd3 cells and ROI perim.
% v2: never ran v1. cell segment with Kathryn alg.
% v3: Lx. hi mag 20x. no overlays. allow cut tissue. xy coord for perim of tissue, ROI, and cd3 in zone Lx, L0, L1, L2.
% v4: Lx. overlays. no cut. xy coord for perim of ROI, and cd3 separate necro in zone Lx, L0, L1, L2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assigns working directory WDIR according to VARARGIN
for k=1:nargin
 	switch varargin{k}
        case '-inputdir'
            if(nargin>=(k+1)) %makes sure there is an input after "-inputdir"
                wdir = varargin{k+1}; %assigns input after "-inputdir" to WDIR
            end
 		case '-flibyer'
 		otherwise
 	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% slide params

path = wdir;    %image
[~, savename] = fileparts(wdir);  % create file identifier
splitsep = '<:\|:>';

apiloc = strcat('/gne/research/apps/gslideviewer/',getenv('GSLIDEMODE'),'/projectonlayer');
capture = 'capture.bash';

%% my params
resdir = ('/gne/home/hungj4/proj/nanozoomer/JZ_10166_cd8count');  % parent dir for all results
mkdir(resdir);
cd(resdir);

savefile2 = strcat(savename,'_stats_Lx.xls');  %slide summary stats

%maskname = 'L0_L1_L2_fix';				%overlay name in gSlide viewer
maskname = 'L0_L1_L2_var';
maskname1 = 'cells_w_necro';
maskname2 = 'cd3';

%clear overlay in gSlide if it exists
opts = struct();
opts.REMOVE = 1;
opts.LAYERNAME = maskname;
OverlayAPI(path,opts);

opts.LAYERNAME = maskname1;
OverlayAPI(path,opts);

opts.LAYERNAME = maskname2;
OverlayAPI(path,opts);

opts.LAYERNAME = 'L0_L1_L2';
OverlayAPI(path,opts);

UID = 50;
tmpfile = sprintf('/tmp/eraseme%s.jpg', UID);
lowmagpower = 2.5;  % magnification power for alignment
himagpower = 20;    % magnification power for texture analysis, normally 10x
scanmag = 20;		% normally 20x.

% control settings
control_w1 = '232'; % this is autotuned later
control_w2 = '10';	% normally 5
control_x1 = '200';
control_x2 = 'd6';
control_x3 = '2000';
control_x4 = 'd5';
expt_id = '1';

% initialize running totals of areas and cell counts
L0_tissue_area = double(0);
L1_tissue_area = double(0);
L2_tissue_area = double(0);
L0_cell_count = double(0);			L0_cell_area = double(0);
L1_cell_count = double(0);			L1_cell_area = double(0);
L2_cell_count = double(0);			L2_cell_area = double(0);
L0_cd3_cell_count = double(0);		L0_cd3_cell_area = double(0);
L1_cd3_cell_count = double(0);		L1_cd3_cell_area = double(0);
L2_cd3_cell_count = double(0);		L2_cd3_cell_area = double(0);
L0_tissue_area_no_necro = double(0);
L1_tissue_area_no_necro = double(0);
L2_tissue_area_no_necro = double(0);
L0_cell_count_no_necro = double(0);			L0_cell_area_no_necro = double(0);
L1_cell_count_no_necro = double(0);			L1_cell_area_no_necro = double(0);
L2_cell_count_no_necro = double(0);			L2_cell_area_no_necro = double(0);
L0_cd3_cell_count_no_necro = double(0);		L0_cd3_cell_area_no_necro = double(0);
L1_cd3_cell_count_no_necro = double(0);		L1_cd3_cell_area_no_necro = double(0);
L2_cd3_cell_count_no_necro = double(0);		L2_cd3_cell_area_no_necro = double(0);

Lx_tissue_area = double(0);
Lx_cell_count = double(0);			Lx_cell_area = double(0);
Lx_cd3_cell_count = double(0);		Lx_cd3_cell_area = double(0);
Lx_trichrome_area = double(0);

Lx_tissue_area_no_necro = double(0);
Lx_cell_count_no_necro = double(0);			Lx_cell_area_no_necro = double(0);
Lx_cd3_cell_count_no_necro = double(0);		Lx_cd3_cell_area_no_necro = double(0);
Lx_trichrome_area_no_necro = double(0);

% create default structuring element
d300 = strel('disk', 150);
d200 = strel('disk', 100);
d150 = strel('disk', 75);
d100 = strel('disk', 50);
d50 = strel('disk', 25);
d40 = strel('disk', 20);
d30 = strel('disk', 15);
d26 = strel('disk', 13);
d22 = strel('disk', 11);
d20 = strel('disk', 10);
d14 = strel('disk', 7);
d10 = strel('disk', 5);	% old d9
d8 = strel('disk', 4);	% old d7
%sq4 = strel('square', 4);
%d3 = strel('square', 3);
d2 = strel('square', 2);

d7 =  logical([ 0 0 1 1 1 0 0;
    0 1 1 1 1 1 0;
    1 1 1 1 1 1 1;
    1 1 1 1 1 1 1;
    1 1 1 1 1 1 1;
    0 1 1 1 1 1 0;
    0 0 1 1 1 0 0]);

d6 =  logical([ 0 1 1 1 1 0;
    1 1 1 1 1 1;
    1 1 1 1 1 1;
    1 1 1 1 1 1;
    1 1 1 1 1 1;
    0 1 1 1 1 0]);

d5 = logical([ 0 1 1 1 0;
    1 1 1 1 1;
    1 1 1 1 1;
    1 1 1 1 1;
    0 1 1 1 0]);

d4 = logical([ 0 1 1 0;
    1 1 1 1;
    1 1 1 1;
    0 1 1 0 ]);
    
d3 = logical([ 1 1 1;
    1 1 1 ;
    1 1 1 ]);

	display('start - get low-res images');
	tCP1 = tic;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get low-res images

%% Slide info, parse JSON data format into Matlab struct, ASTRUCT
opts = struct();
opts.INFO = 1;
opts.MAGPOWER = lowmagpower;
[status, astruct] = CaptureAPI(path,opts);

%% get slide annots - assumes all channels have same info
slideannots = getSlideAnnot(path);

%% get CHAN1 image
tmpfile = getImage(path, lowmagpower);
ss_bigmont = imread(tmpfile);

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Get low-res images takes %d seconds.\n\n', timeCP1);

%% convert background (& necrosis) to white
backgrd = 233;		% can find by looking at image with imtool
control_w2x = '10';	% more open better at this step

g1 = rgb2gray(ss_bigmont);
g2 = uint8(255 * mat2gray(g1));	

lplev1 = (backgrd - str2num(control_w2x)) / 255;
lplev2 = (backgrd + str2num(control_w2x)) / 255;
g3a = xor(im2bw(g2, lplev1), im2bw(g2, lplev2));

xz1 = bwareaopen(g3a, eval(control_x1));	 	% 100, 1000, 2000
xz1 = imclose(xz1, eval(control_x2));  			% d14, d20
xz1 = bwareaopen(xz1, eval(control_x3));
%  	xz1 = imclose(xz1, eval(control_x4)); 
xz1 = ~(bwareaopen(~xz1, 10000));	% closes medium holes

temp_bigmont = ss_bigmont;

	display('start - get ROIs');
	tCP1 = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create mask of ROIS and remove ones named "del", keep ones called "pos"
%%   add mask to capture butterfly ROI

AllRegions = VectorAPI(path, 'GETALL'); %get all ROI ID's to be stamped

del_rois = false(size(ss_bigmont(:,:,1)));  	%blank image to be labeled at 2.5x
add_rois = false(size(ss_bigmont(:,:,1)));  	%blank image to be labeled at 2.5x
fly_rois = false(size(ss_bigmont(:,:,1)));  	%blank image to be labeled at 2.5x
necro_rois = false(size(ss_bigmont(:,:,1)));	%blank image to be labeled at 2.5x
cut_rois = false(size(ss_bigmont(:,:,1)));		%blank image to be labeled at 2.5x

for j = 1 : (size(AllRegions, 2))
        opts = struct();
        opts.MAGPOWER = lowmagpower;
        opts.REGIONID = AllRegions(j).ID;    
        opts.MASK = 1;

        [status, bigmont_mask] = CaptureAPI(path,opts); %get ROI as fragment of whole image at LOWMAGPOWER
		bigmont_mask = imclose(~im2bw(bigmont_mask), d3) ;

        opts.INFO = 1;
        [status, sinfo] = CaptureAPI(path, opts); %get full mag coordinates  of upper left corner of image fragment

        %scale 
        ttop = round(sinfo.t / (scanmag / lowmagpower));  tbottom = ttop + round(size(bigmont_mask, 1)) - 1;  
        tleft = round(sinfo.l / (scanmag / lowmagpower));  tright = tleft + round(size(bigmont_mask, 2)) - 1;

        if ( strcmpi(AllRegions(j).NAME, 'del') || strcmpi(AllRegions(j).NAME, 'exclude') )
            del_rois(ttop : tbottom , tleft : tright) = ...
                imadd( imfill(im2bw(bigmont_mask), 'holes') , del_rois(ttop : tbottom , tleft : tright) );
        elseif strcmpi(AllRegions(j).NAME, 'pos') 
            add_rois(ttop : tbottom , tleft : tright) = ...
                imadd( imfill(im2bw(bigmont_mask), 'holes') , add_rois(ttop : tbottom , tleft : tright) );
        elseif strcmpi(AllRegions(j).NAME, 'necrosis') 
            necro_rois(ttop : tbottom , tleft : tright) = ...
                imadd( imfill(im2bw(bigmont_mask), 'holes') , necro_rois(ttop : tbottom , tleft : tright) );
        elseif strcmpi(AllRegions(j).NAME, 'cut') 
			cut_rois(ttop : tbottom , tleft : tright) = ...
                imadd( imfill(im2bw(bigmont_mask), 'holes') , cut_rois(ttop : tbottom , tleft : tright) );
		else
            fly_rois(ttop : tbottom , tleft : tright) = ...
                imadd( imfill(im2bw(bigmont_mask), 'holes') , fly_rois(ttop : tbottom , tleft : tright) );
        end
end

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Get ROIs takes %d seconds.\n\n', timeCP1);

	display('start - process low-res images');
	tCP1 = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process low-res images, find tissue areas to analyze

% complete void mask
void_mask = xz1 & fly_rois;

temp_bigmont(repmat(~void_mask, [1,1,3])) = 0;		% sets non-void to black
tlevel = 230 / 255;		%  
void_mask = im2bw(temp_bigmont, tlevel);
void_mask = bwareaopen(void_mask, 20);
void_mask = imclose(void_mask, d5);
void_mask = bwareaopen(void_mask, 100);

% grow ROI
L0_mask = fly_rois;

% adjust iterations to get 42 px (lo mag)/337 px(hi mag) [150 um] dilation for L0
L0_mask = imdilate(L0_mask, d40);
L0_mask = imdilate(L0_mask, d22);
L0_mask = imdilate(L0_mask, d22);

% make Lmin mask
Lmin_mask = imdilate(fly_rois, d30);
Lmin_mask = imdilate(Lmin_mask, d26);
Lmin_ring = Lmin_mask & ~fly_rois;


% calc L1 and L2 masks (allows tissue to be cut without false L1 areas)
L1_mask = false(size(ss_bigmont(:,:,1)));  	%blank image to be labeled at 2.5x
L2_mask = false(size(ss_bigmont(:,:,1)));  	%blank image to be labeled at 2.5x

fly_lab = bwlabel(bwareaopen(fly_rois, 500));	% ignore tiny junk
stats = regionprops(fly_lab, 'BoundingBox');

fprintf('The number of ROIs is %d.\n\n', size(stats, 1));

for mask_count = 1 :size(stats, 1)

	mag_ratio = himagpower / lowmagpower;

	% ------  get subset of image from gslide  ------ %
  	start_mask = ismember(fly_lab, mask_count);	% make mask with ONLY one region 	
	start_mask = imfill(start_mask, 'holes');
		temp_mask = bwareaopen(start_mask, (10000/mag_ratio));	% clean up.  MSR10166 cd8 33IW

	% ------  make L1 mask  ------ %
	% L1 width varies depending on area

	tstats = regionprops(temp_mask, 'Area', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity'); 
	a = tstats.MajorAxisLength / 2 ;
	fprintf('Half of major axis length for subset %d is %d.\n\n', mask_count, a);

	% variable L1 width determined by tissue area/length of major axis
	if (a > round(7576/mag_ratio))
		temp_mask = imerode(temp_mask, d200);			% 169 px = 600 um at 2.5x mag
		temp_mask = imerode(temp_mask, d40);
		temp_mask = imerode(temp_mask, d20);
		for ctr = 1:3
			temp_mask = imerode(temp_mask, d26);
		end
		fprintf('L1 width is 600 um.\n\n');
	elseif (a > round(6199/mag_ratio))
		temp_mask = imerode(temp_mask, d200);			% 141 px = 500 um at 2.5x mag
		temp_mask = imerode(temp_mask, d40);
		temp_mask = imerode(temp_mask, d20);
		temp_mask = imerode(temp_mask, d22);	
		fprintf('L1 width is 500 um.\n\n');		
	elseif (a > round(4822/mag_ratio))
		temp_mask = imerode(temp_mask, d200);			% 113 px = 400 um at 2.5x mag
		temp_mask = imerode(temp_mask, d26);
		fprintf('L1 width is 400 um.\n\n');
	elseif (a > round(3444/mag_ratio))
		temp_mask = imerode(temp_mask, d100);			% 84 px = 300 um at 2.5x mag
		temp_mask = imerode(temp_mask, d40);
		temp_mask = imerode(temp_mask, d14);
		temp_mask = imerode(temp_mask, d14);
		fprintf('L1 width is 300 um.\n\n');
	else
		temp_mask = imerode(temp_mask, d40);			% 56 px = 200 um at 2.5x mag
		temp_mask = imerode(temp_mask, d20);
		temp_mask = imerode(temp_mask, d26);
		temp_mask = imerode(temp_mask, d26);
		fprintf('L1 width is 200 um.\n\n');
	end

	L1_ring = start_mask & ~temp_mask;
	L1_mask = L1_mask | L1_ring;

	% ------  make L2 mask  ------ %
	L2_mask = L2_mask | temp_mask;

end  % L1-L2 loop


% find tissue areas to analyze
tlevel = ((backgrd-5) /255);		% was -2 keeps too much backgrd.  -10 removes necrosis
tempmont_mask = ~im2bw(ss_bigmont, tlevel);
	
% density filter to remove non-tissue areas
ss_mask3 = tempmont_mask;
ss_mask3 = imopen(ss_mask3, d2);
ss_mask3 = imclose(ss_mask3, d3); 

ss_mask3 = imclose(ss_mask3, d10);
ss_mask3 = imopen(ss_mask3, d10);				

ss_mask3 = ~(bwareaopen(~ss_mask3, 1000)); 		% closes small holes in tissue
ss_mask3 = bwareaopen(ss_mask3, 1000); 			% removes small pieces of tissue


% fine tune borders of L1 and L2
L1_mask = L1_mask & ss_mask3 ;
L2_mask = L2_mask & ss_mask3 ;


% keep tissue in L0
cutout_mask3 = ss_mask3 & L0_mask;

% remove border tissue < 100 um thickness
cutout_mask3 = cutout_mask3 & ~Lmin_ring;
cutout_mask3 = imclose(cutout_mask3, d100);

% remove voids, but leave necrosis
%necro_mask = necro_mask & ~void_mask;
necro_mask = necro_rois & ~void_mask;
ss_mask3 = ss_mask3 & ~void_mask;
ss_mask3 = bwareaopen(ss_mask3, 10000); % removes pieces of tissue

% leave only area around ROI
ss_mask = ss_mask3 & ~del_rois;
%ss_mask = ss_mask & ~cut_rois;					% cut large tissue to avoid OUT OF MEM
ss_mask = bwareaopen(ss_mask, 100000); 			%  was 20000.  removes large pieces of tissue

% ver Lx. throw away tissue not connected to ROI.
ss_mask = imreconstruct(fly_rois, ss_mask);

ss_mask = bwlabel(im2bw(ss_mask));
stats = regionprops(ss_mask, 'BoundingBox');

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Process low-res images takes %d seconds.\n\n', timeCP1);

	display('start - resize masks');
	tCP1 = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% resize masks to himag
% SS_MASK controls what portions of the slide are pulled back at himag
opts = struct();
opts.INFO = 1;
opts.MAGPOWER = himagpower;
[status, himag] = CaptureAPI(path,opts);

ss_mask_full = imresize(ss_mask, [(astruct.totYpix/himag.scale) (astruct.totXpix/himag.scale)], 'nearest');  %scale analysis SS_MASK to size of full image
ss_mask_full = uint8(ss_mask_full);  % whole tissue mask

add_rois_full = imresize(add_rois, [(astruct.totYpix/himag.scale) (astruct.totXpix/himag.scale)], 'nearest');  %scale to size of full image
%add_rois_full = uint8(add_rois_full);

user_rois_full = imresize(fly_rois, [(astruct.totYpix/himag.scale) (astruct.totXpix/himag.scale)], 'nearest');  %scale to size of full image
%user_rois_full = uint8(user_rois_full);

% ROI with L0 mask
cutout_full = imresize(cutout_mask3, [(astruct.totYpix/himag.scale) (astruct.totXpix/himag.scale)], 'nearest');  %scale to size of full image

necrosis_full = imresize(necro_mask, [(astruct.totYpix/himag.scale) (astruct.totXpix/himag.scale)], 'nearest');  %scale to size of full image
%necrosis_full = uint8(false(size(ss_mask_full, 1) , size(ss_mask_full, 2) , 1 ));

L1_mask_full = imresize(L1_mask, [(astruct.totYpix/himag.scale) (astruct.totXpix/himag.scale)], 'nearest');  %scale to size of full image
L2_mask_full = imresize(L2_mask, [(astruct.totYpix/himag.scale) (astruct.totXpix/himag.scale)], 'nearest');  %scale to size of full image

L0_mask_full = uint8(false(size(ss_mask_full, 1) , size(ss_mask_full, 2) , 1 ));
cellmask_full = uint8(false(size(ss_mask_full, 1) , size(ss_mask_full, 2) , 1 ));
cd3_cell_full = uint8(false(size(ss_mask_full, 1) , size(ss_mask_full, 2) , 1 ));

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Resize masks takes %d seconds.\n\n', timeCP1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process regions from ss_mask separately to make cell mask
% for each distinct region in ss_mask, pull just that image data at himag

fprintf('The number of subsets is %d.\n\n', size(stats, 1));

for mask_count = 1 :size(stats, 1)

	fprintf('Get subset %d from gSlide.\n\n', mask_count);
	tCP1 = tic;

	mcstr = num2str(mask_count);

	% ------  get subset of image from gslide  ------ %
	%mask_count = 1;
	
  	bigmont_mask = ismember(ss_mask_full, mask_count);  %make mask with ONLY one region 	
    hstats = regionprops(bigmont_mask, 'BoundingBox'); %get the bounding box of that region
    
	%calculate real coordinates for image using BoundingBox
	hstats.t = hstats.BoundingBox(1,2) * himag.scale ;
	hstats.l = hstats.BoundingBox(1,1) * himag.scale ;
	hstats.r = (hstats.BoundingBox(1,1) + hstats.BoundingBox(1,3)) * himag.scale ;
	hstats.b = (hstats.BoundingBox(1,2) + hstats.BoundingBox(1,4)) * himag.scale ;

	% ------  when image is too big  ------ %
    %get left half of the image
    opts = struct();
    opts.LEFT = hstats.l;
    opts.RIGHT = floor((hstats.r+hstats.l) / 2);
    opts.TOP = hstats.t;
    opts.BOTTOM = hstats.b;
    opts.MAGPOWER = himagpower;
    opts.OUTPUTFORMAT = 'png';
    [status, bigmont] = CaptureAPI(path,opts);
    
    %get right half and cat
    opts.LEFT = floor((hstats.r+hstats.l) / 2) +1;
    opts.RIGHT = hstats.r;
    [status, bigmont2] = CaptureAPI(path, opts);
    
    bigmont2 = imresize(bigmont2, [size(bigmont, 1) size(bigmont, 2)]) ;    
    bigmont = cat(2, bigmont, bigmont2);
    clear bigmont2;   

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Get subset %d from gSlide takes %d seconds.\n\n', mask_count, timeCP1);

	fprintf('Set non-tissue area to background for subset %d.\n\n', mask_count);
	tCP1 = tic;

	% ------  set non-tissue area to background  ------ %
	% crop add_rois mask 
	add_rois_mask = imcrop(add_rois_full, hstats.BoundingBox);
	add_rois_mask = imresize(add_rois_mask, [size(bigmont, 1) size(bigmont, 2)]);

	% crop user_rois mask 
	user_rois_mask = imcrop(user_rois_full, hstats.BoundingBox);
	user_rois_mask = imresize(user_rois_mask, [size(bigmont, 1) size(bigmont, 2)]);

	% crop L0 cutout mask 
	cutout_mask = imcrop(cutout_full, hstats.BoundingBox);
	cutout_mask = imresize(cutout_mask, [size(bigmont, 1) size(bigmont, 2)]);

	% crop necrosis mask 
	necro_mask = imcrop(necrosis_full, hstats.BoundingBox);
	necro_mask = imresize(necro_mask, [size(bigmont, 1) size(bigmont, 2)]);

	% crop L1 and L2 masks
	L1_mask = imcrop(L1_mask_full, hstats.BoundingBox);
	L1_mask = imresize(L1_mask, [size(bigmont, 1) size(bigmont, 2)]);
	L2_mask = imcrop(L2_mask_full, hstats.BoundingBox);
	L2_mask = imresize(L2_mask, [size(bigmont, 1) size(bigmont, 2)]);

	% make sure mask image matches the size of the image data
    bigmont_mask = imcrop(bigmont_mask, hstats.BoundingBox);
    bigmont_mask = imresize(bigmont_mask, [size(bigmont, 1) size(bigmont, 2)]);

   	mbigmont = bigmont;
   	mbigmont(~repmat(bigmont_mask, [1,1,3])) = backgrd;  % match background, 0 for black

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Set non-tissue area to background for subset %d takes %d seconds.\n\n', mask_count, timeCP1);

	fprintf('Locate cells for subset %d.\n\n', mask_count);
	tCP1 = tic;

	% ------  cell mask  ------ %
	mb_gray = rgb2gray(mbigmont);
	mbg_compl = imcomplement(mb_gray);
	mbg_temp = imtophat(mbg_compl, d10);
	
 	tlevel = (50.5 / 255);		% 50.5 for mbg_compl, 39.5 for mbg_temp.  marked cells
 	tempmont_mask = im2bw(mbg_compl, tlevel);

	tempmont_mask = bwareaopen(tempmont_mask, 20);

	% ------  Kathryn cell segmentation alg  ------ %
	%ci_compl = imcomplement(cell_img);

% 	% convexity-concavity analysis  (not needed.  bimg => ci_compl)
% 	% pre-process image to isolate nuclear signal
% 	comp_top = imcomplement(imtophat(bimg, strel('disk', 35))); % remove background with tophat with large strel and complement
% 	intmask = (comp_top < 240);
% 	intent2 = comp_top;
% 	intent2(~intmask) = 255;
% 	intentC = imclose(intent2, d14); % suppresses dark details smaller than the strel
% 
% 	% create a binary mask that will be used for watershed
% 	dapi_thresh = 240/255;
% 	bimg_mask = ~im2bw(intentC, dapi_thresh);

	% use marker controlled watershed segmentation 
	% see http://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
	%dist = -bwdist(~bimg_mask); % calculates the distance transform using euclidean distance
	dist = -bwdist(~tempmont_mask); % calculates the distance transform using euclidean distance
	mask = imextendedmin(dist, 2);
	dist2 = imimposemin(dist, mask);
	lines = watershed(dist2) == 0;

	% create a mask of cells with line separations due to watershed
	%dapi_mask = bimg_mask & ~lines;
	dapi_mask = tempmont_mask & ~lines;
	dapi_mask = ~bwareaopen(~dapi_mask, 500); 		% fill small holes
	dapi_mask = bwareaopen(dapi_mask, 200); 		% was 800. remove small objects

	% set a size threshold for identifying candidate nuceli clusters
	dapi_clusters = bwareaopen(dapi_mask, 2500);
	%dapi_clusters = bwareaopen(dapi_mask, 800); 	% was 2500.

	% identify concavity points in each candidate cluster
	[dapi_lab, num] = bwlabel(dapi_clusters);

	%concave_points_mask = false(size(bimg)); % preallocate
	%split_lines = false(size(bimg)); % preallocate
	
	concave_points_mask = false(size(tempmont_mask)); % preallocate
	%split_lines = false(size(tempmont_mask)); % preallocate
	split_lines_whole = uint8(false(size(tempmont_mask, 1) , size(tempmont_mask, 2) , 1 ));
	
	warning('off', 'signal:findpeaks:largeMinPeakHeight'); % turn off a warning when no peaks (concavity points) are found

	% in each candidate cluster identify concavity points (if any) by:
	%
	% 1) Measure euclidean distance between the boundary of the dapi signal
	% and the convex hull
	% 2) Identify local peaks in this list of values. These represent
	% concavity points.
	% 3) Perform splitting of clusters with methods based on # of concavity
	% points. For all splits, if area of any resulting fragment is less than a
	% certain threshold don't split.
	%      a) 1 concavity point: measure circularity. if object is too
	%      circular, don't split. Otherwise, split at concavity point and the
	%      midpoint of the dapi signal boundary (based on concavity point
	%      start position).
	%      b) 2 concavity points: split between these two points. Use
	%      bresenham function (from Matlab Mathworks file exchange) to draw 
	%      the rasterized line.
	%      c) > 2 concavity points: use imfindcircles to identify candidate
	%      nuclei. Use the centers of these circles as markers and perform
	%      marker-controlled watershed to generate watershed split lines.

	fprintf('Number of cluster to process for subset %d is %d.\n\n', mask_count, num);

	for ii = 1 : num
		fprintf('Diagnos: Working on cluster %d.\n\n', ii);		% comment out later
		% Create a mask with just the nucleus of interest. Generate the
		% boundary of the dapi signal (nuclear boundary). Generate the boundary of the convex
		% hull of the dapi signal.
		
		%nuc = ismember(dapi_lab, ii);
			nuc_whole = ismember(dapi_lab, ii);
		
			nwstats = regionprops(nuc_whole, 'BoundingBox'); %get the bounding box of that region
	
			%calculate coordinates of nuc_whole using BoundingBox
			nwstats.t = nwstats.BoundingBox(1,2) * himag.scale ;
			nwstats.l = nwstats.BoundingBox(1,1) * himag.scale ;
			nwstats.r = (nwstats.BoundingBox(1,1) + nwstats.BoundingBox(1,3)) * himag.scale ;
			nwstats.b = (nwstats.BoundingBox(1,2) + nwstats.BoundingBox(1,4)) * himag.scale ;
	
			nuc = imcrop(nuc_whole, nwstats.BoundingBox);
			split_lines = false(size(nuc));
			 
		nuc_boundary = bwboundaries(nuc, 8); 
		nuc_boundary = nuc_boundary{1,1};
		hull = bwconvhull(nuc); 
		hull_boundary = bwboundaries(hull, 8); 
		hull_boundary = hull_boundary{1,1};
		minDist_values = zeros(size(nuc_boundary, 1), 1); %pre-allocate

		% For each point on the nuclear boundary, calculate the euclidean
		% distance to all points on the convex hull. Record the minimum value
		% in a vector (minDist_values)
		for jj = 1: size(nuc_boundary, 1)
			ptx = nuc_boundary(jj, 1);
			pty = nuc_boundary(jj, 2);
			xVect = hull_boundary(:, 1);
			yVect = hull_boundary(:, 2);
			distVect = sqrt((xVect-ptx).^2 + (yVect-pty).^2); % euclidean distance
			minDist = min(distVect);
			minDist_values(jj, 1) = minDist;
		end
 
		% Smooth the minimum distance values with a 3 window median filter
		% (helps prevent spurious local peaks). Identify local peaks using
		% findpeaks with a minimum peak height, minimum distance between two
		% peaks, and minimum peak prominence. These settings are set empirically
		% and help prevent spurious local peaks. It can be helpful to plot the
		% minimum distance values and visually identify proper peaks and set
		% options based on that.

		minDist_values = medfilt1(minDist_values, 3);
		[pks, pk_idx] = findpeaks(minDist_values, 'MinPeakHeight', 2, 'MinPeakDistance', 10, 'MinPeakProminence', 3);

		% Use the peak index (pk_idx) to identify the concave_points in the
		% nuclear boundary. Generate a mask of concave points which is useful for
		% visualization and script development.
		concave_points = nuc_boundary(pk_idx, :);
		concave_points_idx = sub2ind(size(concave_points_mask), concave_points(:, 1), concave_points(:, 2));
		concave_points_mask(concave_points_idx) = 1;
 
		% % Splitting decisions
		% fragment size threshold - if any nuclear fragment smaller than threashold don't split
		frag_size_threshold = 200;		% was 500

		if size(concave_points, 1) == 1 % for 1 concavity point
			% if object is too circular, do not perform split
			% Measure circularity based on 4*area*pi/perimeter^2
			% Circularity of a circle = 1. values range 0-1.
			stats = regionprops(nuc, 'Perimeter', 'Area');
			C = (4*stats.Area*pi)./(stats.Perimeter.^2);
			C_threshold = 0.90;
			if C < C_threshold
				% split at the midpoint of the nuclear boundary (midpoint based on starting
				% position defined by concavity point)
				midpoint = floor(size(nuc_boundary, 1)./2); % half the size of the nuc_boundary list
				midpoint_idx = pk_idx + midpoint; % concave point index + midpoint
				nuc_boundary_cat = [nuc_boundary; nuc_boundary]; % double the list so you will continue at the beginning of the list when you reach the end
				split_point = nuc_boundary_cat(midpoint_idx, :);
				[line_x, line_y] = bresenham(concave_points(1,1), concave_points(1,2),...
				 split_point(1,1), split_point(1,2)); % draw the split line with bresenham function
				split_idx = sub2ind(size(split_lines), line_x, line_y);
				split_lines(split_idx) = 1; % add the new split line to the split_lines mask

				% if the split gives a fragment smaller than the threshold, don't split
				nuc_split = nuc & ~split_lines;
				split_lab = bwlabel(nuc_split, 4); % use 4-connectivity because bresenham line is rasterized and objects connect at diagonals
				split_area = regionprops(split_lab, 'Area');
				if min([split_area.Area]) < frag_size_threshold;
					split_lines(split_idx) = 0;
				end
			else
				%continue;		% causes error
			end

		elseif size(concave_points, 1) == 2 % for 2 concavity points
			% split between the two concavity points. draw rasterized line
			% with bresenham function.
			[line_x, line_y] = bresenham(concave_points(1,1), concave_points(1,2), concave_points(2,1), concave_points(2,2));
			split_idx = sub2ind(size(split_lines), line_x, line_y);
			split_lines(split_idx) = 1; % add new split line to split_lines mask

			% if the split gives a fragment smaller than the threshold, don't split
			nuc_split = nuc & ~split_lines;
			split_lab = bwlabel(nuc_split, 4);
			split_area = regionprops(split_lab, 'Area');
			if min([split_area.Area]) < frag_size_threshold;
				split_lines(split_idx) = 0;
			end
	 
		elseif size(concave_points, 1) > 2 % for > 2 concavity points
			% use imfindcircles to identify candidate nuclei
			[centers, radii] = imfindcircles(nuc, [10 30]); 	%radius range 10-30 (set empirically)
			%[centers1, radii] = imfindcircles(nuc, [5 15]);		%
			%[centers2, radii] = imfindcircles(nuc, [16 30]);	%
			%centers = [ centers1; centers2 ];

			if size(centers, 1) > 0			% skip if no centers found
				centers = round(centers); 	% round the values to get integers
				circle_mask = false(size(nuc));
				centers_idx = sub2ind(size(nuc), centers(:, 2), centers(:, 1));
				circle_mask(centers_idx) = 1;
				circles = imdilate(circle_mask, strel('disk', 5));
				% Use the dilated centers of the circles as markers in
				% marker-controlled watershed. Use the resulting watershed lines
				% as split lines.
				watershed_lines = watershed(~circles) == 0;		% if wshed==0, wshed_lines==T
				watershed_lines = nuc & watershed_lines;
				split_lines = split_lines | watershed_lines; % add watershed lines to split lines mask

				% if the split gives a fragment smaller than the threshold, don't split
				nuc_split = nuc & ~split_lines;
				split_lab = bwlabel(nuc_split, 4);				% 4-conn objs
				split_area = regionprops(split_lab, 'Area');
				if min([split_area.Area]) < frag_size_threshold;
					split_lines = split_lines & ~watershed_lines;
				end
			end

		end	% splitting if
 
 		% ------  add split_lines fragment to split_lines image  ------ %
		sl_top = floor(nwstats.t / himag.scale);  sl_bottom = sl_top + size(split_lines, 1) - 1 ;  
		sl_left = floor(nwstats.l / himag.scale);  sl_right = sl_left + size(split_lines, 2) - 1 ;
		split_lines8 = uint8(split_lines);
		try
			split_lines_whole(sl_top : sl_bottom , sl_left : sl_right) = imadd(split_lines8, split_lines_whole(sl_top : sl_bottom , sl_left : sl_right));	
		catch
			% skip splitting cluster when "Subscript indices must either be real positive integers or logicals."
			% don't want to bother debugging this.
			fprintf('Error while working on cluster %d. Skip it.\n\n', ii);		% comment out later
		end
 
	end		% for loop

	% dilate split lines so that you can separate with 8-connectivity (no diagonals connected)
	%split_lines = imdilate(split_lines, true(2));

	tlev = 0.5 / 255;			% all non-zero labels are put in mask
	split_lines_f = im2bw(split_lines_whole, tlev);	
	split_lines_f = imdilate(split_lines_f, true(2));
 
	%cell_mask = dapi_mask & ~split_lines;
	cell_mask = dapi_mask & ~split_lines_f;
	cell_mask = imopen(cell_mask, d2);
	cell_mask = bwareaopen(cell_mask, 40);			% was 20
	cell_mask = imfill(cell_mask, 'holes');
	
	junk_mask = bwareaopen(cell_mask, 100000);	% was 100k.  reject artifacts
	cell_mask = cell_mask & ~junk_mask;

	clear dapi_mask split_lines_f junk_mask;

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Locate cells for subset %d takes %d seconds.\n\n', mask_count, timeCP1);

	fprintf('Locate cd3 cells for subset %d.\n\n', mask_count);
	tCP1 = tic;

	% ------  cd3 cells  ------ %
	red_mask = uint8(mbigmont(:,:,1));
	blu_mask = uint8(mbigmont(:,:,3));
	tlevel = 150 / 255;			% 
	not_blu_mask = ~im2bw(blu_mask, tlevel);
	not_blu_mask = ~bwareaopen(~not_blu_mask, 200); % fill small holes
	not_blu_mask= bwareaopen(not_blu_mask, 50);  	% remove specks

	% separates marked cells from regular cells.  more accurate marked area
	ring_mask = imdilate(bwperim(not_blu_mask), d2);
		%cell_mask_save = cell_mask;
	cell_mask = cell_mask & ~ring_mask ;
	cell_mask = bwareaopen(cell_mask, 50);  	% remove specks
	candidate_cells = imreconstruct(not_blu_mask, cell_mask);	

	dab_lab = bwlabel(candidate_cells);
	cd3_blu_stats = regionprops(dab_lab, blu_mask , 'MeanIntensity' );
    cd3_blu_list = [cd3_blu_stats.MeanIntensity];
		cd3_red_stats = regionprops(dab_lab, red_mask , 'MeanIntensity' );
		cd3_red_list = [cd3_red_stats.MeanIntensity];

	id_cd3 = find((cd3_blu_list - cd3_red_list) < 30) ;    	% select not blue cells
    cd3_mask = ismember(dab_lab, id_cd3);
    
    more_cells = candidate_cells & ~cd3_mask;
	cell_mask = cell_mask | more_cells;						% blue cells go back
	
	clear blu_mask not_blu_mask ring_mask red_mask candidate_cells more_cells;
	
	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Locate cd3 cells for subset %d takes %d seconds.\n\n', mask_count, timeCP1);

	fprintf('Make Lx, L0, L1, L2 masks for subset %d.\n\n', mask_count);
	tCP1 = tic;

	% ------  make Lx mask  ------ %
	Lx_mask = bigmont_mask & ~cutout_mask;

	% ------  make L0 mask  ------ %
	% L0_mask = bigmont_mask & ~user_rois_mask;
	%L0_mask = new_mask & ~user_rois_mask;
	L0_mask = cutout_mask & ~user_rois_mask;
	L0_mask = bwareaopen(L0_mask, 10000);			% clean up

	% ------  make L1 mask  ------ %
	%L1_mask = user_rois_mask & new_mask & L1_mask;
	L1_mask = user_rois_mask & L1_mask;

	% ------  make L2 mask  ------ %
	%L2_mask = user_rois_mask & new_mask & L2_mask;
	L2_mask = user_rois_mask & L2_mask;

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Make Lx masks for subset %d takes %d seconds.\n\n', mask_count, timeCP1);
	
	fprintf('Calculate stats for subset %d.\n\n', mask_count);
	tCP1 = tic;

	% ------  make Lx dist list  ------ %
	Lx_dab_mask = cd3_mask;
   	Lx_dab_mask(~repmat(Lx_mask, [1,1])) = 0;  	% erase non-Lx labels
   	Lx_dab_mask = imreconstruct(Lx_dab_mask, cd3_mask);
		Lx_necro_mask = Lx_dab_mask;
		Lx_necro_mask(~repmat(necro_mask, [1,1])) = 0;	% erase non-necro cells
		Lx_necro_mask = imreconstruct(Lx_necro_mask, Lx_dab_mask);
		Lx_dab_lab_necro = bwlabel(Lx_necro_mask);
		Lx_dab_necro_stats = regionprops(Lx_dab_lab_necro, mb_gray , 'MeanIntensity' , 'Centroid');
			Lx_not_necro_mask = Lx_dab_mask & ~Lx_necro_mask;	% erase necro cells
			Lx_dab_lab_not_necro = bwlabel(Lx_not_necro_mask);
			Lx_dab_not_necro_stats = regionprops(Lx_dab_lab_not_necro, mb_gray , 'MeanIntensity' , 'Centroid');
	leftover_mask = cd3_mask & ~Lx_dab_mask;

    % ------  calculate Lx cell stats  ------ %
	Lx_tissue_area_subset = sum(sum(Lx_mask));
	Lx_cell_mask = cell_mask & Lx_mask;
	Lx_cell_mask = imreconstruct(Lx_cell_mask, cell_mask);
	Lx_cell_area_subset = sum(sum(Lx_cell_mask));
	[~, Lx_cell_count_subset] = bwlabel(Lx_cell_mask);
		Lx_tissue_area_subset_no_necro = sum(sum(Lx_mask & ~necro_mask));
		Lx_cell_mask_no_necro = Lx_cell_mask & ~necro_mask;
		Lx_cell_mask_no_necro = imreconstruct(Lx_cell_mask_no_necro, Lx_cell_mask);
		Lx_cell_area_subset_no_necro = sum(sum(Lx_cell_mask_no_necro));
		[~, Lx_cell_count_subset_no_necro] = bwlabel(Lx_cell_mask_no_necro);
	Lx_cd3_cell_area_subset = sum(sum(Lx_dab_mask));
	[~, Lx_cd3_cell_count_subset] = bwlabel(Lx_dab_mask);
		Lx_cd3_cell_area_subset_no_necro = sum(sum(Lx_not_necro_mask));
		[~, Lx_cd3_cell_count_subset_no_necro] = bwlabel(Lx_not_necro_mask);

	clear Lx_dab_mask Lx_not_necro_mask Lx_dab_lab_not_necro Lx_necro_mask Lx_dab_lab_necro;

	% ------  make L0 dist list  ------ %
	L0_dab_mask = leftover_mask;
   	L0_dab_mask(~repmat(L0_mask, [1,1])) = 0;  	% erase non-L0 labels
   	L0_dab_mask = imreconstruct(L0_dab_mask, leftover_mask);
		L0_necro_mask = L0_dab_mask;
		L0_necro_mask(~repmat(necro_mask, [1,1])) = 0;	% erase non-necro cells
		L0_necro_mask = imreconstruct(L0_necro_mask, L0_dab_mask);
		L0_dab_lab_necro = bwlabel(L0_necro_mask);
		L0_dab_necro_stats = regionprops(L0_dab_lab_necro, mb_gray , 'MeanIntensity' , 'Centroid');
			L0_not_necro_mask = L0_dab_mask & ~L0_necro_mask;	% erase necro cells
			L0_dab_lab_not_necro = bwlabel(L0_not_necro_mask);
			L0_dab_not_necro_stats = regionprops(L0_dab_lab_not_necro, mb_gray , 'MeanIntensity' , 'Centroid');
	leftover_mask = leftover_mask & ~L0_dab_mask;		

    % ------  calculate L0 stats  ------ %
	L0_tissue_area_subset = sum(sum(L0_mask));
	L0_cell_mask = cell_mask & L0_mask;
	L0_cell_mask = imreconstruct(L0_cell_mask, cell_mask);
	L0_cell_area_subset = sum(sum(L0_cell_mask));
	[~, L0_cell_count_subset] = bwlabel(L0_cell_mask);
		L0_tissue_area_subset_no_necro = sum(sum(L0_mask & ~necro_mask));
		L0_cell_mask_no_necro = L0_cell_mask & ~necro_mask;
		L0_cell_mask_no_necro = imreconstruct(L0_cell_mask_no_necro, L0_cell_mask);
		L0_cell_area_subset_no_necro = sum(sum(L0_cell_mask_no_necro));
		[~, L0_cell_count_subset_no_necro] = bwlabel(L0_cell_mask_no_necro);
	L0_cd3_cell_area_subset = sum(sum(L0_dab_mask));
	[~, L0_cd3_cell_count_subset] = bwlabel(L0_dab_mask);
		L0_cd3_cell_area_subset_no_necro = sum(sum(L0_not_necro_mask));
		[~, L0_cd3_cell_count_subset_no_necro] = bwlabel(L0_not_necro_mask);

	clear L0_dab_mask L0_not_necro_mask L0_dab_lab_not_necro L0_necro_mask L0_dab_lab_necro;

	% ------  make L1 dist list  ------ %
	L1_dab_mask = leftover_mask;
   	L1_dab_mask(~repmat(L1_mask, [1,1])) = 0;  	% erase non-L1 labels
   	L1_dab_mask = imreconstruct(L1_dab_mask, leftover_mask);
		L1_necro_mask = L1_dab_mask;
		L1_necro_mask(~repmat(necro_mask, [1,1])) = 0;	% erase non-necro cells
		L1_necro_mask = imreconstruct(L1_necro_mask, L1_dab_mask);
		L1_dab_lab_necro = bwlabel(L1_necro_mask);
		L1_dab_necro_stats = regionprops(L1_dab_lab_necro, mb_gray , 'MeanIntensity' , 'Centroid');
			L1_not_necro_mask = L1_dab_mask & ~L1_necro_mask;	% erase necro cells
			L1_dab_lab_not_necro = bwlabel(L1_not_necro_mask);
			L1_dab_not_necro_stats = regionprops(L1_dab_lab_not_necro, mb_gray , 'MeanIntensity' , 'Centroid');
	leftover_mask = leftover_mask & ~L1_dab_mask;				

    % ------  calculate L1 stats  ------ %
	L1_tissue_area_subset = sum(sum(L1_mask));
	L1_cell_mask = cell_mask & L1_mask;
	L1_cell_mask = imreconstruct(L1_cell_mask, cell_mask);
	L1_cell_area_subset = sum(sum(L1_cell_mask));
	[~, L1_cell_count_subset] = bwlabel(L1_cell_mask);
		L1_tissue_area_subset_no_necro = sum(sum(L1_mask & ~necro_mask));
		L1_cell_mask_no_necro = L1_cell_mask & ~necro_mask;
		L1_cell_mask_no_necro = imreconstruct(L1_cell_mask_no_necro, L1_cell_mask);
		L1_cell_area_subset_no_necro = sum(sum(L1_cell_mask_no_necro));
		[~, L1_cell_count_subset_no_necro] = bwlabel(L1_cell_mask_no_necro);
	L1_cd3_cell_area_subset = sum(sum(L1_dab_mask));
	[~, L1_cd3_cell_count_subset] = bwlabel(L1_dab_mask);
		L1_cd3_cell_area_subset_no_necro = sum(sum(L1_not_necro_mask));
		[~, L1_cd3_cell_count_subset_no_necro] = bwlabel(L1_not_necro_mask);

	clear L1_dab_mask L1_not_necro_mask L1_dab_lab_not_necro L1_necro_mask L1_dab_lab_necro;

	% ------  make L2 dist list  ------ %
	L2_dab_mask = leftover_mask;
   	L2_dab_mask(~repmat(L2_mask, [1,1])) = 0;  	% erase non-L2 labels
   	L2_dab_mask = imreconstruct(L2_dab_mask, leftover_mask);
		L2_necro_mask = L2_dab_mask;
		L2_necro_mask(~repmat(necro_mask, [1,1])) = 0;	% erase non-necro cells
		L2_necro_mask = imreconstruct(L2_necro_mask, L2_dab_mask);
		L2_dab_lab_necro = bwlabel(L2_necro_mask);
		L2_dab_necro_stats = regionprops(L2_dab_lab_necro, mb_gray , 'MeanIntensity' , 'Centroid');
			L2_not_necro_mask = L2_dab_mask & ~L2_necro_mask;	% erase necro cells
			L2_dab_lab_not_necro = bwlabel(L2_not_necro_mask);
			L2_dab_not_necro_stats = regionprops(L2_dab_lab_not_necro, mb_gray , 'MeanIntensity' , 'Centroid');

    % ------  calculate L2 stats  ------ %
	L2_tissue_area_subset = sum(sum(L2_mask));
	L2_cell_mask = cell_mask & L2_mask;
	L2_cell_mask = imreconstruct(L2_cell_mask, cell_mask);
	L2_cell_area_subset = sum(sum(L2_cell_mask));
	[~, L2_cell_count_subset] = bwlabel(L2_cell_mask);
		L2_tissue_area_subset_no_necro = sum(sum(L2_mask & ~necro_mask));
		L2_cell_mask_no_necro = L2_cell_mask & ~necro_mask;
		L2_cell_mask_no_necro = imreconstruct(L2_cell_mask_no_necro, L2_cell_mask);
		L2_cell_area_subset_no_necro = sum(sum(L2_cell_mask_no_necro));
		[~, L2_cell_count_subset_no_necro] = bwlabel(L2_cell_mask_no_necro);
	L2_cd3_cell_area_subset = sum(sum(L2_dab_mask));
	[~, L2_cd3_cell_count_subset] = bwlabel(L2_dab_mask);
		L2_cd3_cell_area_subset_no_necro = sum(sum(L2_not_necro_mask));
		[~, L2_cd3_cell_count_subset_no_necro] = bwlabel(L2_not_necro_mask);

	clear L2_dab_mask L2_dab_mask L2_not_necro_mask L2_dab_lab_not_necro L2_necro_mask L2_dab_lab_necro;
	clear leftover_mask; 


	%% API call to get image info, including microns per pixel calibration
	opts.INFO = 1;    
	[status, sinfo] = CaptureAPI(path, opts);
	AS = (sinfo.umperpixelAtMag) ^ 2 ;   % sq. microns of each pixel
	fprintf('The value of AS is %d .\n\n', AS);

	fprintf('The Lx tissue area for subset %d is %d .\n\n', mask_count, Lx_tissue_area_subset*AS);
	fprintf('The Lx cell area for subset %d is %d .\n\n', mask_count, Lx_cell_area_subset*AS);
	fprintf('The Lx cell count for subset %d is %d .\n\n', mask_count, Lx_cell_count_subset);
	fprintf('The Lx cd3 cell area for subset %d is %d .\n\n', mask_count, Lx_cd3_cell_area_subset*AS);
	fprintf('The Lx cd3 cell count for subset %d is %d .\n\n', mask_count, Lx_cd3_cell_count_subset);
	fprintf('The Lx tissue area without necrosis for subset %d is %d .\n\n', mask_count, Lx_tissue_area_subset_no_necro*AS);
	fprintf('The Lx cell area without necrosis for subset %d is %d .\n\n', mask_count, Lx_cell_area_subset_no_necro*AS);
	fprintf('The Lx cell count without necrosis for subset %d is %d .\n\n', mask_count, Lx_cell_count_subset_no_necro);
	fprintf('The Lx cd3 cell area without necrosis for subset %d is %d .\n\n', mask_count, Lx_cd3_cell_area_subset_no_necro*AS);
	fprintf('The Lx cd3 cell count without necrosis for subset %d is %d .\n\n', mask_count, Lx_cd3_cell_count_subset_no_necro);	

	fprintf('The L0 tissue area for subset %d is %d .\n\n', mask_count, L0_tissue_area_subset*AS);
	fprintf('The L0 cell area for subset %d is %d .\n\n', mask_count, L0_cell_area_subset*AS);
	fprintf('The L0 cell count for subset %d is %d .\n\n', mask_count, L0_cell_count_subset);
	fprintf('The L0 cd3 cell area for subset %d is %d .\n\n', mask_count, L0_cd3_cell_area_subset*AS);
	fprintf('The L0 cd3 cell count for subset %d is %d .\n\n', mask_count, L0_cd3_cell_count_subset);
	fprintf('The L0 tissue area without necrosis for subset %d is %d .\n\n', mask_count, L0_tissue_area_subset_no_necro*AS);
	fprintf('The L0 cell area without necrosis for subset %d is %d .\n\n', mask_count, L0_cell_area_subset_no_necro*AS);
	fprintf('The L0 cell count without necrosis for subset %d is %d .\n\n', mask_count, L0_cell_count_subset_no_necro);
	fprintf('The L0 cd3 cell area without necrosis for subset %d is %d .\n\n', mask_count, L0_cd3_cell_area_subset_no_necro*AS);
	fprintf('The L0 cd3 cell count without necrosis for subset %d is %d .\n\n', mask_count, L0_cd3_cell_count_subset_no_necro);	

	fprintf('The L1 tissue area for subset %d is %d .\n\n', mask_count, L1_tissue_area_subset*AS);
	fprintf('The L1 cell area for subset %d is %d .\n\n', mask_count, L1_cell_area_subset*AS);
	fprintf('The L1 cell count for subset %d is %d .\n\n', mask_count, L1_cell_count_subset);
	fprintf('The L1 cd3 cell area for subset %d is %d .\n\n', mask_count, L1_cd3_cell_area_subset*AS);
	fprintf('The L1 cd3 cell count for subset %d is %d .\n\n', mask_count, L1_cd3_cell_count_subset);
	fprintf('The L1 tissue area without necrosis for subset %d is %d .\n\n', mask_count, L1_tissue_area_subset_no_necro*AS);
	fprintf('The L1 cell area without necrosis for subset %d is %d .\n\n', mask_count, L1_cell_area_subset_no_necro*AS);
	fprintf('The L1 cell count without necrosis for subset %d is %d .\n\n', mask_count, L1_cell_count_subset_no_necro);
	fprintf('The L1 cd3 cell area without necrosis for subset %d is %d .\n\n', mask_count, L1_cd3_cell_area_subset_no_necro*AS);
	fprintf('The L1 cd3 cell count without necrosis for subset %d is %d .\n\n', mask_count, L1_cd3_cell_count_subset_no_necro);	

	fprintf('The L2 tissue area for subset %d is %d .\n\n', mask_count, L2_tissue_area_subset*AS);
	fprintf('The L2 cell area for subset %d is %d .\n\n', mask_count, L2_cell_area_subset*AS);
	fprintf('The L2 cell count for subset %d is %d .\n\n', mask_count, L2_cell_count_subset);
	fprintf('The L2 cd3 cell area for subset %d is %d .\n\n', mask_count, L2_cd3_cell_area_subset*AS);
	fprintf('The L2 cd3 cell count for subset %d is %d .\n\n', mask_count, L2_cd3_cell_count_subset);
	fprintf('The L2 tissue area without necrosis for subset %d is %d .\n\n', mask_count, L2_tissue_area_subset_no_necro*AS);
	fprintf('The L2 cell area without necrosis for subset %d is %d .\n\n', mask_count, L2_cell_area_subset_no_necro*AS);
	fprintf('The L2 cell count without necrosis for subset %d is %d .\n\n', mask_count, L2_cell_count_subset_no_necro);
	fprintf('The L2 cd3 cell area without necrosis for subset %d is %d .\n\n', mask_count, L2_cd3_cell_area_subset_no_necro*AS);
	fprintf('The L2 cd3 cell count without necrosis for subset %d is %d .\n\n', mask_count, L2_cd3_cell_count_subset_no_necro);	

	% ------  subset header  ------ %
	header = sprintf('\n');
	header = sprintf('%s^Lx tissue area, sq. um.', header);
	header = sprintf('%s^Lx cell area, sq. um.^Lx cells, count', header);
	header = sprintf('%s^Lx cd3 cell area, sq. um.^Lx cd3 cells, count', header);
	header = sprintf('%s^L0 tissue area, sq. um.', header);
	header = sprintf('%s^L0 cell area, sq. um.^L0 cells, count', header);
	header = sprintf('%s^L0 cd3 cell area, sq. um.^L0 cd3 cells, count', header);
	header = sprintf('%s^L1 tissue area, sq. um.', header);
	header = sprintf('%s^L1 cell area, sq. um.^L1 cells, count', header);
	header = sprintf('%s^L1 cd3 cell area, sq. um.^L1 cd3 cells, count', header);
	header = sprintf('%s^L2 tissue area, sq. um.', header);
	header = sprintf('%s^L2 cell area, sq. um.^L2 cells, count', header);
	header = sprintf('%s^L2 cd3 cell area, sq. um.^L2 cd3 cells, count', header);

	header = sprintf('%s^Lx tissue area no necro, sq. um.', header);
	header = sprintf('%s^Lx cell area no necro, sq. um.^Lx cells no necro, count', header);
	header = sprintf('%s^Lx cd3 cell area no necro, sq. um.^Lx cd3 cells no necro, count', header);
	header = sprintf('%s^L0 tissue area no necro, sq. um.', header);
	header = sprintf('%s^L0 cell area no necro, sq. um.^L0 cells no necro, count', header);
	header = sprintf('%s^L0 cd3 cell area no necro, sq. um.^L0 cd3 cells no necro, count', header);
	header = sprintf('%s^L1 tissue area no necro, sq. um.', header);
	header = sprintf('%s^L1 cell area no necro, sq. um.^L1 cells no necro, count', header);
	header = sprintf('%s^L1 cd3 cell area no necro, sq. um.^L1 cd3 cells no necro, count', header);
	header = sprintf('%s^L2 tissue area no necro, sq. um.', header);
	header = sprintf('%s^L2 cell area no necro, sq. um.^L2 cells no necro, count', header);
	header = sprintf('%s^L2 cd3 cell area no necro, sq. um.^L2 cd3 cells no necro, count', header);

	logger = sprintf('\n');     
	logger = sprintf('%s^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f', logger ,...
		Lx_tissue_area_subset * AS ,...
		Lx_cell_area_subset * AS , Lx_cell_count_subset ,...
		Lx_cd3_cell_area_subset * AS , Lx_cd3_cell_count_subset ,...
		L0_tissue_area_subset * AS ,...
		L0_cell_area_subset * AS , L0_cell_count_subset ,...
		L0_cd3_cell_area_subset * AS , L0_cd3_cell_count_subset ,...
		L1_tissue_area_subset * AS ,...
		L1_cell_area_subset * AS , L1_cell_count_subset ,...
		L1_cd3_cell_area_subset * AS , L1_cd3_cell_count_subset ,...
		L2_tissue_area_subset * AS ,...
		L2_cell_area_subset * AS , L2_cell_count_subset ,...
		L2_cd3_cell_area_subset * AS , L2_cd3_cell_count_subset );

	logger = sprintf('%s^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f', logger ,...
		Lx_tissue_area_subset_no_necro * AS ,...
		Lx_cell_area_subset_no_necro * AS , Lx_cell_count_subset_no_necro ,...
		Lx_cd3_cell_area_subset_no_necro * AS , Lx_cd3_cell_count_subset_no_necro ,...
		L0_tissue_area_subset_no_necro * AS ,...
		L0_cell_area_subset_no_necro * AS , L0_cell_count_subset_no_necro ,...
		L0_cd3_cell_area_subset_no_necro * AS , L0_cd3_cell_count_subset_no_necro ,...
		L1_tissue_area_subset_no_necro * AS ,...
		L1_cell_area_subset_no_necro * AS , L1_cell_count_subset_no_necro ,...
		L1_cd3_cell_area_subset_no_necro * AS , L1_cd3_cell_count_subset_no_necro ,...
		L2_tissue_area_subset_no_necro * AS ,...
		L2_cell_area_subset_no_necro * AS , L2_cell_count_subset_no_necro ,...
		L2_cd3_cell_area_subset_no_necro * AS , L2_cd3_cell_count_subset_no_necro );

	% ------  add newline at the end  ------ %
	header = sprintf('%s\n', header);
	logger = sprintf('%s\n', logger);     

	% ------  write files  ------ %
	cd(resdir)
	savefile_sub = strcat(savename,'_sub',num2str(mask_count),'_stats.xls');  % subset summary stats
	[cvf] = fopen(savefile_sub,'a+');
	fwrite(cvf, header);
	fwrite(cvf, logger);
	fclose(cvf);

	Lx_tissue_area = Lx_tissue_area_subset + Lx_tissue_area;
	Lx_cell_area = Lx_cell_area_subset + Lx_cell_area;
	Lx_cell_count = Lx_cell_count_subset + Lx_cell_count;
	Lx_cd3_cell_area = Lx_cd3_cell_area_subset + Lx_cd3_cell_area;
	Lx_cd3_cell_count = Lx_cd3_cell_count_subset + Lx_cd3_cell_count;
		Lx_tissue_area_no_necro = Lx_tissue_area_subset_no_necro + Lx_tissue_area_no_necro;
		Lx_cell_area_no_necro = Lx_cell_area_subset_no_necro + Lx_cell_area_no_necro;
		Lx_cell_count_no_necro = Lx_cell_count_subset_no_necro + Lx_cell_count_no_necro;
		Lx_cd3_cell_area_no_necro = Lx_cd3_cell_area_subset_no_necro + Lx_cd3_cell_area_no_necro;
		Lx_cd3_cell_count_no_necro = Lx_cd3_cell_count_subset_no_necro + Lx_cd3_cell_count_no_necro;	

	L0_tissue_area = L0_tissue_area_subset + L0_tissue_area;
	L0_cell_area = L0_cell_area_subset + L0_cell_area;
	L0_cell_count = L0_cell_count_subset + L0_cell_count;
	L0_cd3_cell_area = L0_cd3_cell_area_subset + L0_cd3_cell_area;
	L0_cd3_cell_count = L0_cd3_cell_count_subset + L0_cd3_cell_count;
		L0_tissue_area_no_necro = L0_tissue_area_subset_no_necro + L0_tissue_area_no_necro;
		L0_cell_area_no_necro = L0_cell_area_subset_no_necro + L0_cell_area_no_necro;
		L0_cell_count_no_necro = L0_cell_count_subset_no_necro + L0_cell_count_no_necro;
		L0_cd3_cell_area_no_necro = L0_cd3_cell_area_subset_no_necro + L0_cd3_cell_area_no_necro;
		L0_cd3_cell_count_no_necro = L0_cd3_cell_count_subset_no_necro + L0_cd3_cell_count_no_necro;	

	L1_tissue_area = L1_tissue_area_subset + L1_tissue_area;
	L1_cell_area = L1_cell_area_subset + L1_cell_area;
	L1_cell_count = L1_cell_count_subset + L1_cell_count;
	L1_cd3_cell_area = L1_cd3_cell_area_subset + L1_cd3_cell_area;
	L1_cd3_cell_count = L1_cd3_cell_count_subset + L1_cd3_cell_count;
		L1_tissue_area_no_necro = L1_tissue_area_subset_no_necro + L1_tissue_area_no_necro;
		L1_cell_area_no_necro = L1_cell_area_subset_no_necro + L1_cell_area_no_necro;
		L1_cell_count_no_necro = L1_cell_count_subset_no_necro + L1_cell_count_no_necro;
		L1_cd3_cell_area_no_necro = L1_cd3_cell_area_subset_no_necro + L1_cd3_cell_area_no_necro;
		L1_cd3_cell_count_no_necro = L1_cd3_cell_count_subset_no_necro + L1_cd3_cell_count_no_necro;	

	L2_tissue_area = L2_tissue_area_subset + L2_tissue_area;
	L2_cell_area = L2_cell_area_subset + L2_cell_area;
	L2_cell_count = L2_cell_count_subset + L2_cell_count;
	L2_cd3_cell_area = L2_cd3_cell_area_subset + L2_cd3_cell_area;
	L2_cd3_cell_count = L2_cd3_cell_count_subset + L2_cd3_cell_count;
		L2_tissue_area_no_necro = L2_tissue_area_subset_no_necro + L2_tissue_area_no_necro;
		L2_cell_area_no_necro = L2_cell_area_subset_no_necro + L2_cell_area_no_necro;
		L2_cell_count_no_necro = L2_cell_count_subset_no_necro + L2_cell_count_no_necro;
		L2_cd3_cell_area_no_necro = L2_cd3_cell_area_subset_no_necro + L2_cd3_cell_area_no_necro;
		L2_cd3_cell_count_no_necro = L2_cd3_cell_count_subset_no_necro + L2_cd3_cell_count_no_necro;	

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Calculate stats for subset %d takes %d seconds.\n\n', mask_count, timeCP1);

	fprintf('Add subset %d to main image.\n\n', mask_count);
	tCP1 = tic;

	% ------  add grayscale fragment to full image  ------ %
	ttop = round(hstats.t / himag.scale);  tbottom = ttop + size(mbigmont, 1) - 1;  
    tleft = round(hstats.l / himag.scale);  tright = tleft + size(mbigmont, 2) - 1;

	L0_mask8 = uint8(L0_mask);
	L0_mask_full(ttop : tbottom , tleft : tright ) = imadd(L0_mask8, L0_mask_full(ttop : tbottom , tleft : tright ));	

	cell_mask8 = uint8(cell_mask);
	cellmask_full(ttop : tbottom , tleft : tright ) = imadd(cell_mask8, cellmask_full(ttop : tbottom , tleft : tright ));	

	cd3_mask8 = uint8(cd3_mask);
	cd3_cell_full(ttop : tbottom , tleft : tright ) = imadd(cd3_mask8, cd3_cell_full(ttop : tbottom , tleft : tright ));	

	clear L0_mask8 L1_mask8 L2_mask8 cell_mask8 cd3_mask8 necro_mask8;

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Add subset %d to main image takes %d seconds.\n\n', mask_count, timeCP1);

	fprintf('Write out XY coord of cd3 cells in Lx, L0, L1, L2 for subset %d.\n\n', mask_count);
	tCP1 = tic;

    ttop = round(hstats.t / himag.scale); 	tbottom = ttop + round(size(mbigmont, 1)) - 1;  
    tleft = round(hstats.l / himag.scale);  tright = tleft + round(size(mbigmont, 2)) - 1;

	% ------  write out XY coord of cd3 cells  ------ %
	shift_vector = [tleft, ttop] ;
	Lx_dab_necro_xy = cat(1, Lx_dab_necro_stats.Centroid) ;
	L0_dab_necro_xy = cat(1, L0_dab_necro_stats.Centroid) ;
	L1_dab_necro_xy = cat(1, L1_dab_necro_stats.Centroid) ;
	L2_dab_necro_xy = cat(1, L2_dab_necro_stats.Centroid) ;

	Lx_dab_not_necro_xy = cat(1, Lx_dab_not_necro_stats.Centroid) ;
	L0_dab_not_necro_xy = cat(1, L0_dab_not_necro_stats.Centroid) ;
	L1_dab_not_necro_xy = cat(1, L1_dab_not_necro_stats.Centroid) ;
	L2_dab_not_necro_xy = cat(1, L2_dab_not_necro_stats.Centroid) ;

	nothing = [NaN, NaN];
	
	if (size(Lx_dab_necro_xy, 1) > 0)
		Lx_dab_necro_xy = bsxfun(@plus, Lx_dab_necro_xy, shift_vector) ;
	else
		Lx_dab_necro_xy = nothing;
	end
	if (size(L0_dab_necro_xy, 1) > 0)
		L0_dab_necro_xy = bsxfun(@plus, L0_dab_necro_xy, shift_vector) ;
	else
		L0_dab_necro_xy = nothing;
	end
	if (size(L1_dab_necro_xy, 1) > 0)
		L1_dab_necro_xy = bsxfun(@plus, L1_dab_necro_xy, shift_vector) ;
	else
		L1_dab_necro_xy = nothing;
	end
	if (size(L2_dab_necro_xy, 1) > 0)
		L2_dab_necro_xy = bsxfun(@plus, L2_dab_necro_xy, shift_vector) ;
	else
		L2_dab_necro_xy = nothing;
	end

	if (size(Lx_dab_not_necro_xy, 1) > 0)
		Lx_dab_not_necro_xy = bsxfun(@plus, Lx_dab_not_necro_xy, shift_vector) ;
	else
		Lx_dab_not_necro_xy = nothing;
	end
	if (size(L0_dab_not_necro_xy, 1) > 0)
		L0_dab_not_necro_xy = bsxfun(@plus, L0_dab_not_necro_xy, shift_vector) ;
	else
		L0_dab_not_necro_xy = nothing;
	end
	if (size(L1_dab_not_necro_xy, 1) > 0)
		L1_dab_not_necro_xy = bsxfun(@plus, L1_dab_not_necro_xy, shift_vector) ;
	else
		L1_dab_not_necro_xy = nothing;
	end
	if (size(L2_dab_not_necro_xy, 1) > 0)
		L2_dab_not_necro_xy = bsxfun(@plus, L2_dab_not_necro_xy, shift_vector) ;
	else
		L2_dab_not_necro_xy = nothing;
	end

	savefile_Lx_necro_xy = strcat(savename,'_sub',num2str(mask_count),'_Lx_necro_xy.csv');  % subset cd3 xy coord	
	savefile_L0_necro_xy = strcat(savename,'_sub',num2str(mask_count),'_L0_necro_xy.csv');  % subset cd3 xy coord
	savefile_L1_necro_xy = strcat(savename,'_sub',num2str(mask_count),'_L1_necro_xy.csv');  % subset cd3 xy coord
	savefile_L2_necro_xy = strcat(savename,'_sub',num2str(mask_count),'_L2_necro_xy.csv');  % subset cd3 xy coord

	savefile_Lx_not_necro_xy = strcat(savename,'_sub',num2str(mask_count),'_Lx_not_necro_xy.csv');  % subset cd3 xy coord	
	savefile_L0_not_necro_xy = strcat(savename,'_sub',num2str(mask_count),'_L0_not_necro_xy.csv');  % subset cd3 xy coord
	savefile_L1_not_necro_xy = strcat(savename,'_sub',num2str(mask_count),'_L1_not_necro_xy.csv');  % subset cd3 xy coord
	savefile_L2_not_necro_xy = strcat(savename,'_sub',num2str(mask_count),'_L2_not_necro_xy.csv');  % subset cd3 xy coord

	csvwrite(savefile_Lx_necro_xy, Lx_dab_necro_xy) ;
	csvwrite(savefile_L0_necro_xy, L0_dab_necro_xy) ;
	csvwrite(savefile_L1_necro_xy, L1_dab_necro_xy) ;
	csvwrite(savefile_L2_necro_xy, L2_dab_necro_xy) ;

	csvwrite(savefile_Lx_not_necro_xy, Lx_dab_not_necro_xy) ;
	csvwrite(savefile_L0_not_necro_xy, L0_dab_not_necro_xy) ;
	csvwrite(savefile_L1_not_necro_xy, L1_dab_not_necro_xy) ;
	csvwrite(savefile_L2_not_necro_xy, L2_dab_not_necro_xy) ;

	clear Lx_dab_necro_xy L0_dab_necro_xy L1_dab_necro_xy L2_dab_necro_xy;
	clear Lx_dab_not_necro_xy L0_dab_not_necro_xy L1_dab_not_necro_xy L2_dab_not_necro_xy;

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Write out XY coord of cd3 cells in Lx, L0, L1, L2 for subset %d takes %d seconds.\n\n', mask_count, timeCP1);

	fprintf('Write out YX coord for ROI perimeter for subset %d.\n\n', mask_count);
	tCP1 = tic;

	% ------  write out XY coord of L1/ROI perimeter  ------ %
	flip_vector = [ttop, tleft] ;						% offset when not flipping YX
	%perim_mask = user_rois_mask & new_mask ;
	perim_mask = user_rois_mask & bigmont_mask ;		% MSR10166
	perim_mask = bwareaopen(perim_mask, 200000); 		% remove med to large tissue
	perim_all = bwboundaries(perim_mask, 'noholes') ;
	
	for j = 1 :size(perim_all, 1)
		if (j == 1)
			perim_xy = perim_all{j} ; 
		else
			perim_xy = [perim_xy; perim_all{j}] ;
		end
	end
		
	down_sample_freq = 10 ;					% keeps 1/freq perim points
	perim_xy = perim_xy (1:down_sample_freq:end,:) ;
	perim_xy = bsxfun(@plus, perim_xy, flip_vector) ;
	% plot(perim_xy(:,2), perim_xy(:,1), 'r.') ;
	
	savefile_perim_xy = strcat(savename,'_sub',num2str(mask_count),'_roi_perim_xy.csv');  % subset perim xy coord
	
	csvwrite(savefile_perim_xy, perim_xy) ;

	clear perim_xy;

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Write out YX coord for ROI perimeter for subset %d takes %d seconds.\n\n', mask_count, timeCP1);

end  % process section loop

fprintf('The Lx tissue area is %d .\n\n', Lx_tissue_area*AS);
fprintf('The Lx cell area is %d .\n\n', Lx_cell_area*AS);
fprintf('The Lx cell count is %d .\n\n', Lx_cell_count);
fprintf('The Lx cd3 cell area is %d .\n\n', Lx_cd3_cell_area*AS);
fprintf('The Lx cd3 cell count is %d .\n\n', Lx_cd3_cell_count);
fprintf('The Lx tissue area without necrosis is %d .\n\n', Lx_tissue_area_no_necro*AS);
fprintf('The Lx cell area without necrosis is %d .\n\n', Lx_cell_area_no_necro*AS);
fprintf('The Lx cell count without necrosis is %d .\n\n', Lx_cell_count_no_necro);
fprintf('The Lx cd3 cell area without necrosis is %d .\n\n', Lx_cd3_cell_area_no_necro*AS);
fprintf('The Lx cd3 cell count without necrosis is %d .\n\n', Lx_cd3_cell_count_no_necro);

fprintf('The L0 tissue area is %d .\n\n', L0_tissue_area*AS);
fprintf('The L0 cell area is %d .\n\n', L0_cell_area*AS);
fprintf('The L0 cell count is %d .\n\n', L0_cell_count);
fprintf('The L0 cd3 cell area is %d .\n\n', L0_cd3_cell_area*AS);
fprintf('The L0 cd3 cell count is %d .\n\n', L0_cd3_cell_count);
fprintf('The L0 tissue area without necrosis is %d .\n\n', L0_tissue_area_no_necro*AS);
fprintf('The L0 cell area without necrosis is %d .\n\n', L0_cell_area_no_necro*AS);
fprintf('The L0 cell count without necrosis is %d .\n\n', L0_cell_count_no_necro);
fprintf('The L0 cd3 cell area without necrosis is %d .\n\n', L0_cd3_cell_area_no_necro*AS);
fprintf('The L0 cd3 cell count without necrosis is %d .\n\n', L0_cd3_cell_count_no_necro);

fprintf('The L1 tissue area is %d .\n\n', L1_tissue_area*AS);
fprintf('The L1 cell area is %d .\n\n', L1_cell_area*AS);
fprintf('The L1 cell count is %d .\n\n', L1_cell_count);
fprintf('The L1 cd3 cell area is %d .\n\n', L1_cd3_cell_area*AS);
fprintf('The L1 cd3 cell count is %d .\n\n', L1_cd3_cell_count);
fprintf('The L1 tissue area without necrosis is %d .\n\n', L1_tissue_area_no_necro*AS);
fprintf('The L1 cell area without necrosis is %d .\n\n', L1_cell_area_no_necro*AS);
fprintf('The L1 cell count without necrosis is %d .\n\n', L1_cell_count_no_necro);
fprintf('The L1 cd3 cell area without necrosis is %d .\n\n', L1_cd3_cell_area_no_necro*AS);
fprintf('The L1 cd3 cell count without necrosis is %d .\n\n', L1_cd3_cell_count_no_necro);

fprintf('The L2 tissue area is %d .\n\n', L2_tissue_area*AS);
fprintf('The L2 cell area is %d .\n\n', L2_cell_area*AS);
fprintf('The L2 cell count is %d .\n\n', L2_cell_count);
fprintf('The L2 cd3 cell area is %d .\n\n', L2_cd3_cell_area*AS);
fprintf('The L2 cd3 cell count is %d .\n\n', L2_cd3_cell_count);
fprintf('The L2 tissue area without necrosis is %d .\n\n', L2_tissue_area_no_necro*AS);
fprintf('The L2 cell area without necrosis is %d .\n\n', L2_cell_area_no_necro*AS);
fprintf('The L2 cell count without necrosis is %d .\n\n', L2_cell_count_no_necro);
fprintf('The L2 cd3 cell area without necrosis is %d .\n\n', L2_cd3_cell_area_no_necro*AS);
fprintf('The L2 cd3 cell count without necrosis is %d .\n\n', L2_cd3_cell_count_no_necro);

	fprintf('Write stats file start.\n\n');
	tCP1 = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% log annotations and results
header = sprintf('^^^^^^^');

loc = strcmp(attrib, 'barcode text');
logger = sprintf('%s', char(tvalue(loc)));

logger = sprintf('%s^%s^=hyperlink(RC[-1],"gSlide")^', logger, char(strcat('http://resgslideprd01.gene.com/nano/gslideviewer.cgi?PATH=', path)));

tval = char(tvalue(strcmp(attrib, 'Group')));
if isempty(tval)
    logger = sprintf('%s^^^^', logger); else logger = sprintf('%s%s^^^^', logger, tval); end

tval = char(tvalue(strcmp(attrib, 'Slide Number')));
if ~isempty(tval) 
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Slide Number_TXT', header); 

tval = char(tvalue(strcmp(attrib, 'Block')));
if ~isempty(tval)
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Block_TXT', header); 

tval = char(tvalue(strcmp(attrib, 'Experiment')));
if ~isempty(tval)
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Experiment_TXT', header); 

tval = char(tvalue(strcmp(attrib, 'Tissue Type Name')));
if ~isempty(tval)
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Tissue Type Name_TXT', header); 

tval = char(tvalue(strcmp(attrib, 'Antibody')));
if ~isempty(tval)
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Antibody_TXT', header); 

tval = char(tvalue(strcmp(attrib, 'Group')));
if ~isempty(tval)
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Group_TXT', header); 

tval = char(tvalue(strcmp(attrib, 'Animal')));
if ~isempty(tval)
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Animal_TXT', header); 

tval = char(tvalue(strcmp(attrib, 'Treatment')));
if ~isempty(tval)
    logger = sprintf('%s^%s', logger, tval); else logger = sprintf('%s^', logger); end
header = sprintf('%s^Treatment_TXT', header); 

% ------  add results  ------ %
header = sprintf('%s^Lx tissue area, sq. um.', header);
header = sprintf('%s^Lx cell area, sq. um.^Lx cells, count', header);
header = sprintf('%s^Lx cd3 cell area, sq. um.^Lx cd3 cells, count', header);
header = sprintf('%s^L0 tissue area, sq. um.', header);
header = sprintf('%s^L0 cell area, sq. um.^L0 cells, count', header);
header = sprintf('%s^L0 cd3 cell area, sq. um.^L0 cd3 cells, count', header);
header = sprintf('%s^L1 tissue area, sq. um.', header);
header = sprintf('%s^L1 cell area, sq. um.^L1 cells, count', header);
header = sprintf('%s^L1 cd3 cell area, sq. um.^L1 cd3 cells, count', header);
header = sprintf('%s^L2 tissue area, sq. um.', header);
header = sprintf('%s^L2 cell area, sq. um.^L2 cells, count', header);
header = sprintf('%s^L2 cd3 cell area, sq. um.^L2 cd3 cells, count', header);

header = sprintf('%s^Lx tissue area no necro, sq. um.', header);
header = sprintf('%s^Lx cell area no necro, sq. um.^Lx cells no necro, count', header);
header = sprintf('%s^Lx cd3 cell area no necro, sq. um.^Lx cd3 cells no necro, count', header);
header = sprintf('%s^L0 tissue area no necro, sq. um.', header);
header = sprintf('%s^L0 cell area no necro, sq. um.^L0 cells no necro, count', header);
header = sprintf('%s^L0 cd3 cell area no necro, sq. um.^L0 cd3 cells no necro, count', header);
header = sprintf('%s^L1 tissue area no necro, sq. um.', header);
header = sprintf('%s^L1 cell area no necro, sq. um.^L1 cells no necro, count', header);
header = sprintf('%s^L1 cd3 cell area no necro, sq. um.^L1 cd3 cells no necro, count', header);
header = sprintf('%s^L2 tissue area no necro, sq. um.', header);
header = sprintf('%s^L2 cell area no necro, sq. um.^L2 cells no necro, count', header);
header = sprintf('%s^L2 cd3 cell area no necro, sq. um.^L2 cd3 cells no necro, count', header);

% ------  add newline at the start  ------ %
%logger = sprintf('%s\n', logger);     

logger = sprintf('%s^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f', logger ,...
	Lx_tissue_area * AS ,...
	Lx_cell_area * AS , Lx_cell_count ,...
	Lx_cd3_cell_area * AS , Lx_cd3_cell_count ,...
	L0_tissue_area * AS ,...
	L0_cell_area * AS , L0_cell_count ,...
	L0_cd3_cell_area * AS , L0_cd3_cell_count ,...
	L1_tissue_area * AS ,...
	L1_cell_area * AS , L1_cell_count ,...
	L1_cd3_cell_area * AS , L1_cd3_cell_count ,...
	L2_tissue_area * AS ,...
	L2_cell_area * AS , L2_cell_count ,...
	L2_cd3_cell_area * AS , L2_cd3_cell_count );

logger = sprintf('%s^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f^%f', logger ,...
	Lx_tissue_area_no_necro * AS ,...
	Lx_cell_area_no_necro * AS , Lx_cell_count_no_necro ,...
	Lx_cd3_cell_area_no_necro * AS , Lx_cd3_cell_count_no_necro ,...
	L0_tissue_area_no_necro * AS ,...
	L0_cell_area_no_necro * AS , L0_cell_count_no_necro ,...
	L0_cd3_cell_area_no_necro * AS , L0_cd3_cell_count_no_necro ,...
	L1_tissue_area_no_necro * AS ,...
	L1_cell_area_no_necro * AS , L1_cell_count_no_necro ,...
	L1_cd3_cell_area_no_necro * AS , L1_cd3_cell_count_no_necro ,...
	L2_tissue_area_no_necro * AS ,...
	L2_cell_area_no_necro * AS , L2_cell_count_no_necro ,...
	L2_cd3_cell_area_no_necro * AS , L2_cd3_cell_count_no_necro );

% ------  add newline at the end  ------ %
header = sprintf('%s\n', header);
logger = sprintf('%s\n', logger);     

% ------  write files  ------ %
cd(resdir)
[cvf] = fopen(savefile2,'a+');
fwrite(cvf, header);
fwrite(cvf, logger);
fclose(cvf);

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Write stats file takes %d seconds.\n\n', timeCP1);

	fprintf('Write overlays to gSlide start.\n\n');
	tCP1 = tic;

%% convert lp mask to bw
tlev = 0.5 / 255;			% all non-zero labels are put in mask
hr_temp_mask_f = im2bw(L0_mask_full, tlev);
hr_perim_mask_f = imdilate(bwperim(hr_temp_mask_f),strel('disk',5));

hr_temp_mask_f = im2bw(L1_mask_full, tlev);
hr_perim_mask_f = hr_perim_mask_f | imdilate(bwperim(hr_temp_mask_f),strel('disk',5));

hr_temp_mask_f = im2bw(L2_mask_full, tlev);
hr_perim_mask_f = hr_perim_mask_f | imdilate(bwperim(hr_temp_mask_f),strel('disk',5));

hr_cell_mask_f =  im2bw(cellmask_full, tlev);
hr_cd3_cell_f =  im2bw(cd3_cell_full, tlev);
hr_necrosis_f =  im2bw(necrosis_full, tlev);

	% stop out of memory errors
	clear L0_mask_full L1_mask_full L2_mask_full ;
	clear cellmask_full ;
	clear cd3_cell_full ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add masks to overlay image in gSlide

%%%%%APPARENT SIZE LIMITATION in RGB2IND, stick to 10x or less

% tissue outline mask - red

toverlay = cat(3, uint8(hr_perim_mask_f).* 255, uint8(hr_perim_mask_f) .* 0, uint8(hr_perim_mask_f).* 0);
toverlay = imresize(toverlay, .5, 'nearest');
alpha = uint8(rgb2gray(toverlay));
alpha(alpha>0) = 200;

opts = struct();
opts.LAYERNAME = maskname;
opts.PROJECT = 1;
opts.SRCIMAGE = toverlay;
opts.SRCMAG = himagpower * .5;
opts.TOP = 0;
opts.LEFT = 0;
opts.BOTTOM = astruct.totYpix;
opts.RIGHT = astruct.totXpix;
opts.ALPHA = alpha;

[stat, res] = OverlayAPI(path,opts);

% cell_mask - blue,  necrosis = orange
mc_overlay = cat(3, uint8(hr_necrosis_f).* 255, uint8(hr_necrosis_f) .* 128, uint8(hr_cell_mask_f).* 255);
mc_overlay = imresize(mc_overlay, .5, 'nearest');
alpha = uint8(rgb2gray(mc_overlay));
alpha(alpha>0) = 200;

opts = struct();
opts.LAYERNAME = maskname1;
opts.PROJECT = 1;
opts.SRCIMAGE = mc_overlay;
opts.SRCMAG = himagpower * .5;
opts.TOP = 0;
opts.LEFT = 0;
opts.BOTTOM = astruct.totYpix;
opts.RIGHT = astruct.totXpix;
opts.ALPHA = alpha;

[stat, res] = OverlayAPI(path,opts);

% cd3 cells - green
mc_overlay = cat(3, uint8(hr_cd3_cell_f).* 0, uint8(hr_cd3_cell_f).* 255, uint8(hr_cd3_cell_f).* 0);
mc_overlay = imresize(mc_overlay, .5, 'nearest');
alpha = uint8(rgb2gray(mc_overlay));
alpha(alpha>0) = 200;

opts = struct();
opts.LAYERNAME = maskname2;
opts.PROJECT = 1;
opts.SRCIMAGE = mc_overlay;
opts.SRCMAG = himagpower * .5;
opts.TOP = 0;
opts.LEFT = 0;
opts.BOTTOM = astruct.totYpix;
opts.RIGHT = astruct.totXpix;
opts.ALPHA = alpha;

[stat, res] = OverlayAPI(path,opts);


% stop out of memory errors
	fun = 0; fun2 = 0;   % avoid load warnings
	clear toverlay;
	clear mc_overlay;

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Write overlays to gSlide takes %d seconds.\n\n', timeCP1);


%---------------------- nested functions ----------------------------------
%
% Function to get gslide viewer info

    function slideannots = getSlideAnnot (path_name)

        % get slide annotations, store attributes in ATTRIB, values in TVALUE
        [s, slideannots] = system( sprintf( ...
			strcat('/gne/research/apps/gslideviewer/',getenv('GSLIDEMODE'),'/projectonlayer/getannots.bash -path "%s"'),...
        path_name));
        slideannots = regexp(slideannots, splitsep, 'split');
        slideannots = reshape(slideannots, [], 1);

        [attrib, tvalue] = strtok(slideannots, ':');
        for i = 1:size(tvalue, 1)
            tvalue(i,1) = cellstr(tvalue{i,1}(2:end));
        end

    end % of getSlideAnnot
%--------------------------------------------------------------------------
% Function to get image from dB

    function returnfile = getImage (path_name, magnification)

    % get dimensions of image
    %     get info from API
	cmd = sprintf('%s/%s -path "%s" -magpower %f -info ',apiloc,capture,path_name, magnification);
	% fprintf('cmd:%s\n',cmd);  % debug from Steve
    [s, slideinfo] = system(cmd);

    Xpix=str2double(regexprep(regexp(slideinfo,'totXpix=\w*','match'),'totXpix=',''));
    Ypix=str2double(regexprep(regexp(slideinfo,'totYpix=\w*','match'),'totYpix=',''));

    % get low-res or hi-res image
	cmd = sprintf('%s/%s -path "%s" -magpower %f -l 0 -b %f -r %f -t 0 -outfile "%s"',...
	        apiloc, capture, path_name, magnification, Ypix, Xpix, tmpfile);
	% fprintf('cmd:%s\n',cmd);  % debug from Steve
    system(cmd);
        
    returnfile = tmpfile;

    end % of getImage
%--------------------------------------------------------------------------


end  % of function JZiai_MSR10166_cd3



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Secondary functions

%--------------------------------------------------------------------------

function [imgout] = rgb_thresh(tempimg, redl, redt, greenl, greent, bluel, bluet)
	% MASKOUT  applies RGB color thresholding accoring to inputs and returns
	% combined mask

	%% separates color channels
	imgred = tempimg(:,:,1);
	imggreen = tempimg(:,:,2);
	imgblue = tempimg(:,:,3);
	clear img

	%% Foci threshold
	%'bot' masks are the lower threshold, 'top' masks are upper.  Union of
	%all three color channels determine mask that
	%is passed to processing step
	if redl > 0
		redbot = ~im2bw(imgred, redl/255);
	else
		redbot = zeros(size(imgred));
	end

	if redt <255
		redtop = ~im2bw(imgred, redt/255);
	else
		redtop = ones(size(imgred));
	end
	redmask = xor (redtop, redbot);
	clear redbot redtop;

	if greenl > 0
		greenbot = ~im2bw(imggreen, greenl/255);
	else
		greenbot = zeros(size(imggreen));
	end

	if greent <255
		greentop = ~im2bw(imggreen, greent/255);
	else
		greentop = ones(size(imggreen));
	end
	greenmask = xor(greentop, greenbot);
	clear greenbot greentop;

	if bluel > 0
		bluebot = ~im2bw(imgblue, bluel/255);
	else
		bluebot = zeros(size(imgblue));
	end

	if bluet <255
		bluetop = ~im2bw(imgblue, bluet/255);
	else
		bluetop = ones(size(imgblue));
	end
	bluemask = xor(bluetop, bluebot);
	clear bluetop bluebot;

	imgout  = redmask & greenmask & bluemask ;
end

%--------------------------------------------------------------------------

function [imgout] = DABT_auto(imgin, thresh)
	% DABT_auto accepts an rgb image as 'imgin', and a numerical value between
	% 0 and 255 as 'thresh', and returns a binary image that corresponds to the
	% pixels whose "brown-ness" exceeds the value given in 'thresh'.  The
	% blue-normalized formula below is taken from Brey et. al. The Journal of 
	% Histochemistry & Cytochemistry, Volume 51(5): 575-584, 2003.

	%separate the individual color channels
	imgred = uint16(imgin(:,:,1));
	imggreen = uint16(imgin(:,:,2));
	imgblue = uint16(imgin(:,:,3)); clear imgin;
 
	denom = imgred + imggreen + imgblue;
	numer = imgblue*255;

	mask = imdivide( numer , denom); clear numer denom;

	%convert mask back to a uint8
	mask2 = uint8(mask);

	%returns a binary mask corresponding to pixels with sufficient brown
	imgout = (mask2 < thresh);
end

%--------------------------------------------------------------------------

function [imgout] = HEMT_auto(imgin, thresh)
	% HEMT_auto accepts an rgb image as 'imgin', and a numerical value between
	% 0 and 255 as 'thresh', and returns a binary image that corresponds to the
	% pixels whose "blue-ness" exceeds the value given in 'thresh'.

	imgred = uint16(imgin(:,:,1));
	imggreen = uint16(imgin(:,:,2));
	imgblue = uint16(imgin(:,:,3));

	averg = (imgred + imggreen) / 2 ;
	mask = imsubtract( imgblue , averg);

	mask2 = uint8(mask);
	imgout = (mask2 > thresh);
end

%--------------------------------------------------------------------------

function [x, y]= bresenham(x1,y1,x2,y2)
	% Matlab optmized version of Bresenham line algorithm. No loops.
	% Format:
	%               [x y]=bham(x1,y1,x2,y2)
	%
	% Input:
	%               (x1,y1): Start position
	%               (x2,y2): End position
	%
	% Output:
	%               x y: the line coordinates from (x1,y1) to (x2,y2)
	%
	% Usage example:
	%               [x y]=bham(1,1, 10,-5);
	%               plot(x,y,'or');
	x1=round(x1); x2=round(x2);
	y1=round(y1); y2=round(y2);
	dx=abs(x2-x1);
	dy=abs(y2-y1);
	steep=abs(dy)>abs(dx);
	if steep t=dx;dx=dy;dy=t; end

	%The main algorithm goes here.
	if dy==0 
		q=zeros(dx+1,1);
	else
		q=[0;diff(mod([floor(dx/2):-dy:-dy*dx+floor(dx/2)]',dx))>=0];
	end

	%and ends here.

	if steep
		if y1<=y2 y=[y1:y2]'; else y=[y1:-1:y2]'; end
		if x1<=x2 x=x1+cumsum(q);else x=x1-cumsum(q); end
	else
		if x1<=x2 x=[x1:x2]'; else x=[x1:-1:x2]'; end
		if y1<=y2 y=y1+cumsum(q);else y=y1-cumsum(q); end
	end
end

%--------------------------------------------------------------------------
