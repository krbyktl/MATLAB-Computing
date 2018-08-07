%kleo1,05/07/2018,R2017b
%Title: An Analysis of U.S. Lung and Bronchus Cancer Rates since 1975
%https://seer.cancer.gov/faststats/selections.php?#Output

%% plan
%Will first plot the sum of mortality and incidence over time to describe
%the general rate of people affected by lung and bronchus cancer. Then, the
%maximum rate will be found from this dataset, and a cut-off of the points
%before and right after this time will be plotted and fitted to compare
%increasing and decreasing speed of cancer rates.
%Contextual information will be provided to possibly explain any
%correlations and reasoning for the data trends.
%% access data
f = fopen('lung_cancer_data.dat');
C = textscan(f,'%s %s %f %f','HeaderLines',5,'Delimiter',',')
fclose(f)

%% organize data
incidence_rate = C{3}(1:41);
mortality_rate = C{3}(42:82);
y = C{2}(1:41);
y1 = regexprep(y,"'",'');
year = str2double(y1);
affected_rate = incidence_rate + mortality_rate
average_incidence_rate = mean(incidence_rate)
average_mortality_rate = mean(mortality_rate)
average_affected_rate = mean(affected_rate)

%% table
data_table=table(year,incidence_rate,mortality_rate,affected_rate,'VariableNames',{'Year','Incidence','Mortality','Affected'})
writetable(data_table)

%% plot of total rate of lung and bronchus cancer affected people
close all
[critical_point,index] = max(affected_rate)
critical_year = year(index)
set(gcf,'color','w')
plot(year,affected_rate,'--b',critical_year,critical_point,'ro','MarkerSize',15,'LineWidth',3)
xlabel('Year','FontSize',16)
ylabel('Affected (# per 100,000)','FontSize',16)
title({'Number of People Affected by', 'Lung and Bronchus Cancer Since 1975'},'FontSize',20)
l=legend('path of yearly trend','maximum point')
set(l,'FontSize',14)
print('total_affected.jpg','-djpeg','-r400')

%% plot before and after critical year
close all
set(gcf,'color','w')
years_before = year(1:index-1);
affected_rate_before = affected_rate(1:index-1);
[f1,gof1] = fit(years_before,affected_rate_before,'poly2')
subplot(1,2,1)
plot(years_before,affected_rate_before,'co')
hold on
plot(years_before,f1(years_before),':g','LineWidth',2)
xlabel('Year','FontSize',16)
ylabel('Affected (# per 100,000)','FontSize',16)
title({'A. Affected Number','Before Maximum'},'FontSize',14)
l1=legend('raw data','y = -0.0689x^2+275.3x-2.748e+05')
set(l1,'Location','South','FontSize',8)
gof1.sse
S_1 = ((gof1.sse).^(1/2))/(length(years_before)-1)

years_after = year(index+1:41);
affected_rate_after = affected_rate(index+1:41);
[f2,gof2] = fit(years_after,affected_rate_after,'poly2')
subplot(1,2,2)
plot(years_after,affected_rate_after,'co')
hold on
plot(years_after,f2(years_after),':g','LineWidth',2)
xlabel('Year','FontSize',16)
ylabel('Affected (# per 100,000)','FontSize',16)
title({'B. Affected Number','After Maximum'},'FontSize',14)
l2=legend('raw data','y = -0.06258x^2+249.3-2.481e+05')
set(l2,'FontSize',8,'Location','South')
gof2.sse
S_2 = ((gof2.sse).^(1/2))/(length(years_after)-1)

print('before_after_fits.jpg','-djpeg','-r400')

%% Rate Analysis
close all
f1x=differentiate(f1,years_before)
plot(years_before,f1x,'-.r','LineWidth',3)
[f11,gof11]=fit(years_before,f1x,'poly1');
hold on
f2x=differentiate(f2,years_after)
plot(years_after,f2x,'--m','LineWidth',3)
[f21,gof21]=fit(years_after,f2x,'poly1');
grid on
set(gcf,'color','w')
xlabel('Year','FontSize',16)
ylabel('\Delta Affected (# per 100,000/year)','FontSize',16)
title('Increase and Decrease in Affected Rate over 40 Years','FontSize',16)
ylim([-4 4])
l3=legend('y = -0.1378x+275.3','y = -0.1252x+249.3')
set(l3,'FontSize',14)
print('rate_analysis.jpg','-djpeg','-r400')
rate_difference = f11.p1 - f21.p1


%% analysis to text
g = fopen('fit_analysis.txt','w')
fprintf(g,'SSE of A:%f\nSSE of B:%f\nS of A:%f\nS of B:%f\nSlope A:%f\nSlope B:%f\nSlope Difference:%f\n',gof1.sse,gof2.sse,S_1,S_2,f11.p1,f21.p1,rate_difference)
fclose(g)
