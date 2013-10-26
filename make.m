function make()
disp('Compiling...');
p = input('Parallel?(Y/N)','s');
if(strcmpi(p,'y')==1)
    if(exist('OCTAVE_VERSION'))
        disp('Compile: P');
        str = '  -D__PARALLEL__';
    else
        disp('Compile: P');
        str = ' -largeArrayDims CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -D__PARALLEL__';
    end
else
    disp('Compile: NP');
    str = ' -largeArrayDims';
end

eval(['mex patchmatch.cpp' str]);
disp('Done.');
end


