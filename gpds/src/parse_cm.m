load coalmine.txt
cm = cumsum(coalmine);
cm = cm*120./cm(end);
cm = cm(1:end-1);
