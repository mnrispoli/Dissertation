hor_cut='000hor_8atoms_cuttinX2p5Y-6p5.raw'
vert_cut = '000vert_left_line_cuttingX3p6Y-7p00.raw'
disorder = '001QPL_low_pxl_defocus_n30_0_of_200.raw'
round_walls = '001round_walls_8sites_shortX3p6Y-7p00.raw'
vert_dw_walls = '002vert_dw_smallbump_flattop18_leftshifX3p6Y-7p00.raw'
final_wall = '004vert_berlin_wall_smoothboX3p6Y-7p00.raw'

path_to_hologram = vert_cut

fin=fopen(path_to_hologram,'r');
I=fread(fin,768*1024,'ubit1');
I=reshape(I,[1024 768]);
figure(1)
imshow(I)

map=[237/256, 28/256, 36/256; 1, 1, 1; ]
colormap(map)