function JZiai_MSR9005_cd3(varargin) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JZIAI_MSR9005_CD3 processes a bright-field image from one slide. 
%
% Filename: JZiai_MSR9005_cd3.m
% Compile (if necessary, next two lines):
%   mcc -m JZiai_MSR9005_cd3 
% Use (within Matlab):  JZiai_MSR9005_cd3
% -inputdir 'NDP_Pcore2/Genentech_Webslide_server/MSRs/20160509_DDunlap_MSR8957/gslide_filename1.ndpi' 
%
% Output files:   
%   stats in .xls file
%
% Dependencies that must be in the same directory as this file.
%   -- none
%
% v1: new code.  zones from MSR8394.
% v2: tune necro
% v3: tune necro again, tune cd3, custom L1
% v4: split stats for each subset.  xy coord for cd3 cells and ROI perim.  Write overlays.
% v5: copy of JZiai_MSR10166_cd3_Lx
% v6: put v4 cell and cd8 find into v5.
% v7: custom fixed L1 from v3.
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
resdir = ('/gne/home/hungj4/proj/nanozoomer/JZ_9005_cd8count');  % parent dir for all results
mkdir(resdir);
cd(resdir);

savefile2 = strcat(savename,'_stats_Lx.xls');  %slide summary stats

maskname = 'L0_L1_L2';				%overlay name in gSlide viewer
maskname1 = 'cells_w_necro';	 	%overlay name in gSlide viewer
maskname2 = 'cd3'; 					%overlay name in gSlide viewer
maskname3 = 'cd3_border'; 

%clear overlay in gSlide if it exists
opts = struct();
opts.REMOVE = 1;
opts.LAYERNAME = maskname;
OverlayAPI(path,opts);

opts.LAYERNAME = maskname1;
OverlayAPI(path,opts);

opts.LAYERNAME = maskname2;
OverlayAPI(path,opts);

opts.LAYERNAME = maskname3;
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
% complete necrosis mask
tlevel = 200 / 255;
tm = ~im2bw(ss_bigmont, tlevel);
tm = tm & ~xz1;
necro_dist = uint8(bwdist(tm));

distlevel = 3 / 255;
necro_mask = im2bw(necro_dist, distlevel);
necro_mask = imclose(necro_mask, d10);
necro_mask = bwareaopen(necro_mask, 200);
necro_mask = imclose(necro_mask, d20);
necro_mask = ~(bwareaopen(~necro_mask, 5000)); 		% closes med holes in tissue
fr = imerode(fly_rois, d50);
necro_mask = necro_mask & fr;
necro_mask = bwareaopen(necro_mask, 200);

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

	% fixed L1 width determined by pathologist
	if (strcmp(savename,'4MTQ-PL16 - 2016-05-10 08.53.13') | strcmp(savename,'4MU2-PL16 - 2016-05-10 09.14.20'))
		temp_mask = imerode(temp_mask, d200);			% 141 px = 500 um at 2.5x mag
		temp_mask = imerode(temp_mask, d40);
		temp_mask = imerode(temp_mask, d20);
		temp_mask = imerode(temp_mask, d22);	
		fprintf('L1 width is 500 um.\n\n');		
	elseif (strcmp(savename,'4MTG-PL16 - 2016-05-09 15.54.54') | strcmp(savename,'4MTL-PL16 - 2016-05-10 08.44.04') | ...
			strcmp(savename,'4MTT-PL16 - 2016-05-10 08.58.50') | strcmp(savename,'4MTU-PL16 - 2016-05-10 09.00.37') | ...
			strcmp(savename,'4MTX-PL16 - 2016-05-10 09.05.49') | strcmp(savename,'4MU0-PL16 - 2016-05-10 09.11.25') | ...
			strcmp(savename,'4MU1-PL16 - 2016-05-10 09.12.57'))
		temp_mask = imerode(temp_mask, d200);			% 113 px = 400 um at 2.5x mag
		temp_mask = imerode(temp_mask, d26);
		fprintf('L1 width is 400 um.\n\n');
	elseif strcmp(savename,'4MTW-PL16 - 2016-05-10 09.04.09')
		temp_mask = imerode(temp_mask, d100);			% 84 px = 300 um at 2.5x mag
		temp_mask = imerode(temp_mask, d40);
		temp_mask = imerode(temp_mask, d14);
		temp_mask = imerode(temp_mask, d14);
		fprintf('L1 width is 300 um.\n\n');
	elseif strcmp(savename,'4MTY-PL16 - 2016-05-10 09.07.21')
		temp_mask = imerode(temp_mask, d40);			% 56 px = 200 um at 2.5x mag
		temp_mask = imerode(temp_mask, d20);
		temp_mask = imerode(temp_mask, d26);
		temp_mask = imerode(temp_mask, d26);
		fprintf('L1 width is 200 um.\n\n');
	else
		temp_mask = imerode(temp_mask, d200);			% 169 px = 600 um at 2.5x mag
		temp_mask = imerode(temp_mask, d40);
		temp_mask = imerode(temp_mask, d20);
		for ctr = 1:3
			temp_mask = imerode(temp_mask, d26);
		end
		fprintf('L1 width is 600 um.\n\n');
	end

	L1_ring = start_mask & ~temp_mask;
	L1_mask = L1_mask | L1_ring;

	% ------  make L2 mask  ------ %
	L2_mask = L2_mask | temp_mask;

end  % L1-L2 loop


% find tissue areas to analyze
tlevel = ((backgrd-2) /255);		% -10 removes necrosis
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
necro_mask = necro_mask & ~void_mask;
%necro_mask = necro_rois & ~void_mask;
ss_mask3 = ss_mask3 & ~void_mask;
ss_mask3 = bwareaopen(ss_mask3, 10000); % removes pieces of tissue

% leave only area around ROI
ss_mask = ss_mask3 & ~del_rois;
%ss_mask = ss_mask & ~cut_rois;					% cut large tissue to avoid OUT OF MEM
ss_mask = bwareaopen(ss_mask, 20000); 			%  was 100000.  removes large pieces of tissue

% ver Lx. throw away tissue not connected to ROI.
ss_mask = imreconstruct(fly_rois, ss_mask);

ss_mask = bwlabel(im2bw(ss_mask));
stats = regionprops(ss_mask, 'BoundingBox');

% save necro mask for cd8/trichrome
save(strcat(savename, '_necro.mat'),'necro_mask');

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
	
	tlevel = (50.5 / 255);		% was 39.5
	tempmont_mask = im2bw(mbg_compl, tlevel);

	tempmont_mask = bwareaopen(tempmont_mask, 20);

	cell_img = mbigmont;
	cell_img(repmat(~tempmont_mask, [1,1,3])) = 255;

    %% Fast radial transform
    fun = @(block_struct) radial_seg_cells(block_struct.data, d5 , [5 8 ] , 2 , -10 , strel('disk', 14));
    lines = blockproc(cell_img, [2000 2000], fun, 'BorderSize', [30 30]);

	temp = ~bwareaopen(lines, 6000);	% was 6000
    lines = lines & temp;				% avoid seg faults
    lines2 = imopen(lines, d3);
    lines2 = bwareaopen(lines2, 20);
    
    marker = lines & lines2;
    recon_F = imreconstruct(marker, lines2);

	tm_large = bwareaopen(tempmont_mask, 150);
	tm_small = tempmont_mask & ~tm_large;
	tm_large = tm_large & recon_F;
	
	cell_mask = tm_large | tm_small;
	cell_mask = imopen(cell_mask, d2);
	cell_mask = bwareaopen(cell_mask, 40);		% was 20
	cell_mask = imfill(cell_mask, 'holes');

	junk_mask = bwareaopen(cell_mask, 100000);	% was 100k.  reject artifacts
	cell_mask = cell_mask & ~junk_mask;

	clear mbg_compl mbg_temp cell_img lines lines2 marker recon_F tm_small tm_large junk_mask;

	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Locate cells for subset %d takes %d seconds.\n\n', mask_count, timeCP1);

	fprintf('Locate cd3 cells for subset %d.\n\n', mask_count);
	tCP1 = tic;

	% ------  cd3 cells  ------ %
	cell_img = mbigmont;
	cell_img(repmat(~cell_mask, [1,1,3])) = 0;		% only evaluate cells

	red_mask = uint8(cell_img(:,:,1));
	blu_mask = uint8(cell_img(:,:,3));

	tlevel = 150 / 255;			% 
	not_blu_mask = ~im2bw(uint8(mbigmont(:,:,3)), tlevel);
	not_blu_mask = bwareaopen(not_blu_mask, 50);  	% remove specks

	% separates marked cells from regular cells.  more accurate marked area
	ring_mask = imdilate(bwperim(not_blu_mask), d2);
	cell_mask = cell_mask & ~ring_mask ;
	cell_mask = bwareaopen(cell_mask, 50);  	% remove specks

	missing_not_blu = not_blu_mask & ~imreconstruct(cell_mask, not_blu_mask) ;
	missing_not_blu = imfill(missing_not_blu, 'holes');

	candidate_cells = cell_mask | missing_not_blu ;

	dab_lab = bwlabel(candidate_cells);
	cd3_blu_stats = regionprops(dab_lab, blu_mask , 'MeanIntensity' );
    cd3_blu_list = [cd3_blu_stats.MeanIntensity];
		cd3_red_stats = regionprops(dab_lab, red_mask , 'MeanIntensity' );
		cd3_red_list = [cd3_red_stats.MeanIntensity];

	id_cd3 = find((cd3_blu_list - cd3_red_list) < 20) ;    	% select not blue cells
    cd3_mask = ismember(dab_lab, id_cd3);
 
 	cd3_lab = bwlabel(cd3_mask);
	cd3_stats = regionprops(cd3_lab, mb_gray , 'MeanIntensity' );
    cd3_list = [cd3_stats.MeanIntensity];
 
 	id_cd3 = find(cd3_list < 170) ;    	% select not blue cells
    cd3_mask = ismember(cd3_lab, id_cd3);
 	
	clear blu_mask not_blu_mask ring_mask red_mask candidate_cells more_cells;
	
	timeCP1 = toc(tCP1);  % in seconds
	fprintf('Locate cd3 cells for subset %d takes %d seconds.\n\n', mask_count, timeCP1);

	fprintf('Make Lx, L0, L1, L2 masks for subset %d.\n\n', mask_count);
	tCP1 = tic;

	% ------  make Lx mask  ------ %
	Lx_mask = bigmont_mask & ~cutout_mask;

	% ------  make L0 mask  ------ %
	L0_mask = cutout_mask & ~user_rois_mask;
	L0_mask = bwareaopen(L0_mask, 10000);			% clean up

	% ------  make L1 mask  ------ %
	L1_mask = user_rois_mask & L1_mask;

	% ------  make L2 mask  ------ %
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
	perim_mask = user_rois_mask & bigmont_mask ;		% MSR9228
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


end  % of function JZiai_MSR9005_cd3



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

%RADIAL_SEG_CELLS.M is a wrapper function to count cells via radial
%symmetry.  After smoothening by opening and closing with reconstruction by
%the kernel SSMOOTH, for each pixel in image IMGIN, a measure of radial 
%symmetry is calculated at the distances specified in RAD_VECT, and
%constrained by the slack variable RAD_CRITERIA.  Background pixels are
%defined as a sink for watershed transformation as being those left after
%dilating the radial symmetry centers by CELL_KERNEL.
%
%Dependant Matlab toolboxes = Image Processing
%Dependant functions = Radial_Sym_Transform
%
%INPUTS:
%IMGIN is the raw, RGB image.
%
%SSMOOTH is the kernel to use in the initial closing and opening w/recon
% EX:   ssmooth = strel('disk', 3) ;     
%   
%RAD_METRIC radial measure metric; the more negative, the more radial an 
%object needs to be to be kept as a radial center.
% EX:   rad_metric = -25 ;
%
%RAD_VECT vector containing distances to measure radial symmetry from 
%center pixel
% EX:   rad_vect = [4 7 10] ;
%
%RAD_CRITERIA radial strictness parameter; 1 = slack 3 = strict
% EX:   rad_criteria = 2 ;
%
%CELL_DIAMETER kernel to define distance away from radial symmetry centers 
%that are highly likely to be background
% EX:   cell_diameter = strel('disk', 14);
%
%OUTPUT:
%SEGOUT is a binary mask with 0 corresponding to cell borders
%
%
%EXAMPLE:
% lines = radial_seg_cells(rgb_img , d3 , [3 6 9] , 2 , -15 , d9) ;
%
%
function segout = radial_seg_cells (imgin , ssmooth , rad_vect , rad_criteria , rad_metric , cell_kernel )
       
    big_gray = rgb2gray(imgin);  
	%big_gray = imgin;			% in this case, input image already gray
    
    %open by recon
    big_recon = imerode(big_gray, ssmooth);
    big_recon = imreconstruct(big_recon, big_gray);
    
    %close by recon
    big_recon_C = imdilate(big_recon, ssmooth);
    big_recon_C = imreconstruct(imcomplement(big_recon_C), imcomplement(big_recon));
    big_recon_C = imcomplement(big_recon_C);
    
    nuc_bound = imclose(big_recon_C, ssmooth);

    % radial symmetry
    S = Radial_Sym_Transform(double(big_gray), rad_vect , rad_criteria );
    
    nmark = S < rad_metric ;
    nmark = bwareaopen(nmark, 5); %removes radial centers smaller than given size
    
    % find background
    bmask = imdilate(nmark, cell_kernel);    
    bskel = bwmorph(~bmask, 'skel', Inf);
    
    %sobel gradient of preprocessed
    sar = single([1 2 1; 0 0 0; -1 -2 -1]);
    
    H = conv2(single(nuc_bound), sar);
    V = conv2(single(nuc_bound), sar');
    sgradient = sqrt(H.^2 + V.^2);
    sgradient = sgradient(2:end-1, 2:end-1);
    
    sgradient = imimposemin(sgradient, (bskel | nmark) );

    lines = watershed(sgradient);
    segout = lines > 0 ;
    
end

%--------------------------------------------------------------------------

% Radial_Sym_Transform - Loy and Zelinski's fast radial feature detector
%
% An implementation of Loy and Zelinski's fast radial feature detector
%
% Usage: S = fastradial(im, radii, alpha);
%
% Arguments:
%            im    - image to be analysed (after using imread for example:im=imread('lenna.jpg');)
%            radii - array of integer radius values to be processed
%                    suggested radii might be [1 3 5]
%            alpha - radial strictness parameter.
%                    1 - slack, accepts features with bilateral symmetry.
%                    2 - a reasonable compromise.
%                    3 - strict, only accepts radial symmetry.
%                        ... and you can go higher
%
% Returns    S      - Symmetry map.  Bright points with high symmetry are
%                     marked with large positive values. Dark points of
%                     high symmetry marked with large negative values.
%
% To localize points use NONMAXSUPPTS on S, -S or abs(S) depending on
% what you are seeking to find.

% Reference:
% Loy, G.  Zelinsky, A.  Fast radial symmetry for detecting points of
% interest.  IEEE PAMI, Vol. 25, No. 8, August 2003. pp 959-973.

function S = Radial_Sym_Transform(im, radii, alpha)
    if any(radii ~= round(radii)) || any(radii < 1)
        error('radii must be integers and > 1')
    end
    
    [rows,cols]=size(im);
    
    % Use the Sobel masks to get gradients in x and y
    gx = [-1 0 1
          -2 0 2
          -1 0 1];
    gy = gx';
    
    imgx = filter2(gx,im);
    imgy = filter2(gy,im);
    mag = sqrt(imgx.^2 + imgy.^2)+eps; % (+eps to avoid division by 0)
    
    % Normalise gradient values so that [imgx imgy] form unit 
    % direction vectors.
    imgx = imgx./mag;   
    imgy = imgy./mag;
    
    S = zeros(rows,cols);  % Symmetry matrix
    
    [x,y] = meshgrid(1:cols, 1:rows);
    
    for n = radii
	M = zeros(rows,cols);  % Magnitude projection image
	O = zeros(rows,cols);  % Orientation projection image

        % Coordinates of 'positively' and 'negatively' affected pixels
        posx = x + round(n*imgx);
        posy = y + round(n*imgy);
        
        negx = x - round(n*imgx);
        negy = y - round(n*imgy);
        
        % Clamp coordinate values to range [1 rows 1 cols]
        posx( find(posx<1) )    = 1;
        posx( find(posx>cols) ) = cols;
        posy( find(posy<1) )    = 1;
        posy( find(posy>rows) ) = rows;
        
        negx( find(negx<1) )    = 1;
        negx( find(negx>cols) ) = cols;
        negy( find(negy<1) )    = 1;
        negy( find(negy>rows) ) = rows;
        

        % Form the orientation and magnitude projection matrices
        for r = 1:rows
            for c = 1:cols
                O(posy(r,c),posx(r,c)) = O(posy(r,c),posx(r,c)) + 1;
                O(negy(r,c),negx(r,c)) = O(negy(r,c),negx(r,c)) - 1;
                
                M(posy(r,c),posx(r,c)) = M(posy(r,c),posx(r,c)) + mag(r,c);
                M(negy(r,c),negx(r,c)) = M(negy(r,c),negx(r,c)) - mag(r,c);
            end
        end
        
        % Clamp Orientation projection matrix values to a maximum of 
        % +/-kappa,  but first set the normalization parameter kappa to the
        % values suggested by Loy and Zelinski
        if n == 1, kappa = 8; else kappa = 9.9; end
        
        O(find(O >  kappa)) =  kappa;  
        O(find(O < -kappa)) = -kappa;  
        
        % Unsmoothed symmetry measure at this radius value
        F = M./kappa .* (abs(O)/kappa).^alpha;
        
        % Generate a Gaussian of size proportional to n to smooth and spread 
        % the symmetry measure.  The Gaussian is also scaled in magnitude
        % by n so that large scales do not lose their relative weighting.
        A = fspecial('gaussian',[n n], 0.25*n) * n;  
        
        S = S + filter2(A,F);
        
    end  % for each radius
    
    S = S/length(radii);  % Average 
end

%--------------------------------------------------------------------------
