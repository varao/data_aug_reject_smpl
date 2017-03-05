function ll = gpds_lik(l, num_acc, offset)

  ll = sum(log(sigmoid(l(1:num_acc), offset)));
  ll = ll + sum(log(1-sigmoid(l(num_acc+1:end), offset)));
