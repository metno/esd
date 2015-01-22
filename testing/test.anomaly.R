library(esd)

print('Get daily stations: NACD')
y1 <- station(src='ecad',cntr='Norway',nmin=70)

print('Derive monthly values')
y2 <- as.monthly(y1)

print('Derive seasonal values')
y3 <- as.4seasons(y1)
y4 <- as.4seasons(y2)

print('Anomalies of single daily station')
plot(anomaly(subset(y1,is=1)))

print('Anomalies of single monthly station')
plot(anomaly(subset(y2,is=1)))

print('Anomalies of single seasonal station')
plot(anomaly(subset(y3,is=1)))

print('Anomalies of single seasonal station')
plot(anomaly(subset(y4,is=1)))

print('Anomalies of all daily stations')
plot(anomaly(y1))

print('Anomalies of all monthly stations')
plot(anomaly(y2))

print('Anomalies of all seasonal stations')
plot(anomaly(y3))

print('Anomalies of all seasonal stations')
plot(anomaly(y4))
