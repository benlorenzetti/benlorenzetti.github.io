clear;
clc;
run('econdata.m'); % Functions for importing and manipulating data

mort30us = s2m_data(import_fred_csv('..\data\MORTGAGE30US.csv'));   % percent apr
mspus = s2m_data(import_fred_csv('..\data\MSPUS.csv'));             % dollars
rent_index = s2m_data(import_fred_csv('..\data\CUUR0000SEHA.csv')); % unitless
recessions = import_fred_csv('..\data\JHDUSRGDPBR.csv');  % 0 or 1

d = 0.2;
apr_range = months_range(mort30us, mort30us);
apr = stretch_align_index(mort30us, apr_range);
mpr = (1 + apr(2,:)./100).^(1/12);
Rmort = ((1-mpr.^(-360)) ./ (1 - 1./mpr)) / (1-d);
Rreal = 180;
Rth(1,:) = apr(1,:);
Rth(2,:) = (Rreal * Rmort)./(Rmort + Rreal);

norm_year = 2015;
norm_month = 12;
pr_range = months_range(mspus, rent_index);
Vmed = stretch_align_index(mspus, pr_range);
Irent = stretch_align_index(rent_index, pr_range);
ind1 = lookup(Rth(1,:), 12*norm_year + norm_month);
ind2 = lookup(Vmed(1,:), 12*norm_year + norm_month);
scale = Rth(2,ind1) / (Vmed(2,ind2)/(Irent(2,ind2)));
pr_index(1,:) = pr_range(1) : pr_range(2);
pr_index(2,:) = scale * Vmed(2,:)./Irent(2,:);

plot(pr_index(1,:)/12, pr_index(2,:), "linewidth", 2, ';Median Sales Price / City Rent Average;');
hold on;
plot(Rth(1,:)/12, Rth(2,:), "linewidth", 2, ';Mortgage Months Ratio R(M);');
plot(recessions(1,:)/12, 140*recessions(2,:)+40, "k:;US Recessions;");
hold off;
titl = title("US Price-to-Rent Ratio and Federal Interest Rate Policy R(M)");
lege = legend("location","northwest");
ylab = ylabel("Months of Rent");
text(1967, 32, "Data from FRED.stlouisfed.org. MSPUS, CUUR0000SEHA, MORTGAGE30US and JHDUSRGDPBR.");
text(1963, 28.5, "R(M) calculated assuming 30-yr, 20%-down, Rreal=180 months. P/R index normalized to R(M) at (Dec. 2015, 97.3)");
axis([1968, 2023,40,140]);
set(lege, "fontsize", 14);
set(titl, "fontsize", 14);
set(ylab, "fontsize", 14);
