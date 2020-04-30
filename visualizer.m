function visualizer(img, cnn, patchsize)
%VISUALIZER - Visualizes the Similarity
%   visualizer(img, cnn, patchsize)
%   	img - input image (default - 'img' from workspace)
%	cnn - nn-field (default - 'cnn' from workspace)
%	patchsize - patchsize (default - 'patchsize' from workspace)
%
%   Author: Sk. Mohammadul Haque
%   Copyright (c) 2013 Sk. Mohammadul Haque
%


ok = 0;
h2 = 0;

if(nargin<3), patchsize = evalin('base','patchsize'); end;
if(nargin<2), cnn = evalin('base','cnn'); end;
if(nargin<1), img = evalin('base','img'); end;
h = figure;
imshow(uint8(img));
sz = size(img);

if(numel(patchsize)==1), patchsize = [patchsize patchsize]; end;

while true
    try
        % try getting
        while(~ok)
            figure(h);
            [y, x] = ginput(1);
            x = fix(x); y = fix(y);
            if(x>0 && x<=(sz(1)-patchsize(1)) && y>0 && y<=(sz(2)-patchsize(2)))
                ok = 1;
            end
        end
        
        hold off;
        figure(h);
        for k = find(h2~=0)
            delete(h2(k,1));
        end
        hold on;
    
        % get indices
        nx = squeeze(full(cnn(x,y,1,:)))+1;
        ny = squeeze(full(cnn(x,y,2,:)))+1;
        nd = squeeze(full(cnn(x,y,3,:)));
        
        % now draw
        figure(h);
        h2 = zeros(length(nd),1);
        [~,ndi] = sort(nd); 
        ndim = max(ndi);
        for k = 1:length(nd)
            color = [(1-(find(ndi==k,1,'first')/ndim)^0.6) 0.15 0.15];
            h2(k,1) = rectangle('Position',[nx(k) ny(k) patchsize(2) patchsize(1)],'LineWidth', 2,'EdgeColor', color);
        end
        
        % reset
        ok = 0;
     catch dummy
         return;
     end
    
end
end
        
    
    
