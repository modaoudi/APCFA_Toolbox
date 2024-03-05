function  ind = near_value(array, value)
   [~, ind] = min(abs(array - value));
end