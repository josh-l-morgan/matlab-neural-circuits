function sp = nexread2006( filename )

% sp = nexread( filename ) reads a text file exported from NeuroExplorer. It returns a structure sp.
% [filename, pathname]=uigetfile;
% cd = pathname;
fid = fopen( filename , 'rt' );

% Read strings between tabs one at a time until something that is not a label is
% encountered. Once one is found, back up.
nchan=1;
channames{nchan}=fscanf(fid,'%s\t',1);
% Check whether the last read matches a known label format.
while isempty( str2num( channames{nchan} ) )
    nchan=nchan+1;
    % Keep the current file position in case the next read is not a label.
    lastpos=ftell(fid);
    channames{nchan}=fscanf(fid,'%s\t',1);
end

% Hack for new Matlab.
lastpos=lastpos+1;

% The last read was not a label. Back up and drop the last read result.
fseek(fid, lastpos, 'bof');
nchan=nchan-1;
channames=channames(1:end-1);

% Now that the labels have been read, read the numerical data.
sp.data=textscan( fid, '%n', 'delimiter', '\t');
%sp.data=textscan( fid, '%n' );
% Get rid of the cell.
sp.data=sp.data{1};

fclose(fid);

% Find the number of data samples for each channel.
nrows = ( length(sp.data) / nchan );

% Reshape the data by channel.
sp.data=reshape( sp.data, nchan, nrows );
% textscan puts NaN for empty entries. Change them to -1 for compatibility
% with my old code.
sp.data( isnan( sp.data ) ) = -1;
% Transpose because I want nrows by nchannels.
sp.data=sp.data';

sp.channels=channames;
sp.nchannels=nchan;

sp.x=-100*ones(1,nchan);
sp.y=-100*ones(1,nchan);
tmp=regexp(sp.channels,'ch_[0-9][0-9]');
electrodeinds=zeros(size(sp.channels));
for cnt = 1:length(electrodeinds)
    electrodeinds(cnt)=~isempty(tmp{cnt});
end
electrodeinds=find(electrodeinds);
electrodenames=cat(1, sp.channels{ electrodeinds } );
sp.x(electrodeinds)=str2num( electrodenames(:,4) )';
sp.y(electrodeinds)=-str2num( electrodenames(:,5) )';
