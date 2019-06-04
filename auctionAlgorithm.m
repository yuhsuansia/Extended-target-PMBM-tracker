function [Customer2Item,Item2Customer]=auctionAlgorithm(RewardMatrix)

epsilon=0.1; %amount of deviation from the optimal reward

[NofCustomers,NofItems]=size(RewardMatrix);
if (NofCustomers>NofItems)
    error('Number of columns must be greater than or equal to the number of rows');
end

Item2Customer=zeros(1,NofItems);
Customer2Item=zeros(1,NofCustomers);

while ~isempty(find(Customer2Item==0, 1)),
    if (NofItems==1) %if there is only one item
        [~, Item2Customer]=max(RewardMatrix);%Assign the item to the best customer
        Customer2Item(Item2Customer)=1;%Assign the corresponding customer to the item
    else
        for i=1:NofCustomers,
            if ~Customer2Item(i),
                [maxval,maxind]=max(RewardMatrix(i,:));%find maximum element value and its index
                RewardMatrix(i,maxind)=min(RewardMatrix(i,:))-1;%make the maximum minimum to find second maximum
                [secondmaxval,secondmaxind]=max(RewardMatrix(i,:));%find the second maximum value and its index
                RewardMatrix(i,maxind)=maxval; %restore the maximum value

                Customer2Item(i)=maxind; %Assign the customer the item
                if Item2Customer(maxind),%if item is already assigned
                    Customer2Item(Item2Customer(maxind))=0;%unassign the corresponding customer
                end
                Item2Customer(maxind)=i; %Assign the item to the customer
                RewardMatrix(:,maxind)=RewardMatrix(:,maxind)-(maxval-secondmaxval+epsilon);%reduce the item's value
            end
        end
    end
end