function serial_data = parallel_to_serial(parallel_data)
%SERIAL_TO_PARALLEL Function takes Serial Datastream and converts into Nsub
%Parallel Datastreams
%   Detailed explanation goes here
% serial_data (1:length(parallel_data(1,:))*N_sub) = 0;
%     for i = 1:length(parallel_data(1,:))
%         x = i * N_sub - N_sub + 1;
%         y = i * N_sub;
%         serial_data(x:y) = parallel_data(:,i); 
%     end
% 
% end

serial_data = reshape(parallel_data,1,[]);