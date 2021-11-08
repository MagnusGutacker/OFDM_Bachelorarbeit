function parallel_data = serial_to_parallel(serial_data , N_sub  , symbol_size)
%SERIAL_TO_PARALLEL Function takes Serial Datastream and converts into Nsub
%Parallel Datastreams
%   Detailed explanation goes here

parallel_data = zeros(N_sub,int32(length(serial_data)/(N_sub*symbol_size)),symbol_size);
% counter = 1;
%     for i = 1:length(serial_data) / (N_sub * symbol_size)       
%         for j = 1 : N_sub
%             parallel_data(j,i,:) = serial_data(counter:counter+symbol_size-1);
%             counter = counter + symbol_size;
%         end
%     end
parallel_data = reshape(serial_data,symbol_size,N_sub,[]);
parallel_data = permute(parallel_data,[2 3 1]);
end

