%PATCHMATCH - Matches patches in an image and gives an NN-field
%   nnf = patchmatch(img_src, img_dst, patchsize, nnk, searchradius, mask, searchstep, includeself, incomplete, threshold, disttype)
%       img_src - input source image (RGB/Grayscale in double datatype)
%       img_dst - input destination image (RGB/Grayscale in double datatype) (default - img_src)
%       patchsize - size of search patch (default - [5 5])
%       nnk - number of nearest neighbours (default - 5)
%       searchradius - patch search radius (default - 10)
%       mask - patch distance weight mask (size - patchsize X patchsize)(default - ones)
%       searchstep - patch search step (default - 2)
%       includeself - 0/1 (default - 1)
%       incomplete - allow nan values in img - 0/1 (default - 0)
%       threshold - threshold for incomplete image input (default - 0)
%       disttype - distance type, 0 - LARK / 1 - l1 / 2 - l2 / 3 - poisson (default - 2) (LARK not completely implemented)
%       nnf - output nearest neighbour field
%
%   Author: Sk. Mohammadul Haque
%

