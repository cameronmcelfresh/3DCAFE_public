function [gridSmooth] = modeFilter(grid,maskSize)
%modeFilter Functiont to perform mode-based smoothing on the image

gridSmooth = grid;

for c = 1:length(grid)
    for r = 1:length(grid)

        %Row start and end
        if r+maskSize>length(grid)
            rEnd=length(grid);
        else
            rEnd=r+maskSize;
        end

        if r-maskSize<1
            rStart=1;
        else
            rStart=r-maskSize;
        end

        %Column start and end
        if c+maskSize>length(grid)
            cEnd=length(grid);
        else
            cEnd=c+maskSize;
        end

        if c-maskSize<1
            cStart=1;
        else
            cStart=c-maskSize;
        end

        points = grid(rStart:rEnd,cStart:cEnd);

        gridSmooth(r,c)=mode(points(:));

    end
end

end