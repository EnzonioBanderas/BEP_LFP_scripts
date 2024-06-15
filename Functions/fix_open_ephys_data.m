function fix_open_ephys_data(inFileName)
% This function might fix corrupted OpenEphys .continuous files
    fileHeaderSize = 1024;
    blockHeaderSize = 12;
    blockDataLength = 1024;
    recordMarker = uint8([0 1 2 3 4 5 6 7 8 255]);
    
    blockDataSize = 2*blockDataLength;
    recordMarkerSize = length(recordMarker);
    blockSize = blockDataSize + blockHeaderSize + recordMarkerSize;

    outFileName = [inFileName '.recons.continuous'];
    fin = fopen(inFileName,'rb');
    if (fin < 0)
        error('Unable to open input file');
    end
    fout = fopen(outFileName,'wb');
    if (fout < 0)
        fclose(fin);
        error('Unable to open output file');
    end
    
    %copy header;
    D=fread(fin,[1 fileHeaderSize],'uint8');
    fwrite(fout,D,'uint8');
    
    readMarker = zeros(1,recordMarkerSize);
    
    blockPos = ftell(fin);
    
    while(~feof(fin))
        while (~isequal(readMarker,recordMarker))
            readMarker(1:end-1) = readMarker(2:end);
            read = fread(fin,1,'uint8');
            if (feof(fin))
                break;
            end
            readMarker(end)= read;
        end
        nextBlockPos = ftell(fin);
        readSize = nextBlockPos - blockPos;
        if (readSize > recordMarkerSize)
            fseek(fin,blockPos,-1);
            if (readSize > blockSize)
                badBlockSize = readSize - blockSize;
                sprintf('found bad block at %X',...
                    blockPos)
                D=fread(fin,[1 badBlockSize],'uint8');
                fwrite(fout,D,'uint8');
                fwrite(fout,uint8(zeros(1,blockSize - recordMarkerSize - badBlockSize)),'uint8');
                fwrite(fout,recordMarker,'uint8');
            end
            D=fread(fin,[1 blockSize],'uint8');
            fwrite(fout,D,'uint8');
        end
        blockPos = nextBlockPos;
        readMarker = zeros(1,recordMarkerSize);
    end
    
    fclose(fin);
    fclose(fout);
end