source("co_laplace_consistency.R")

#Consistency
#u=0.1 v=0.1 omega is varying
res111 = repfun(N=300, n=23400, u=0.5, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res121 = repfun(N=300, n=23400, u=0.5, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3, gam=0)
res131 = repfun(N=300, n=23400, u=0.5, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3, gam=0)

res.omegac1 = rbind(res111$output, res121$output, res131$output)
write.table(res.omegac1, "example results for omega  c1-u05v05.txt")

res211 = repfun(N=300, n=23400, u=1, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res221 = repfun(N=300, n=23400, u=1,  v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res231 = repfun(N=300, n=23400, u=1, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac2 = rbind(res211$output, res221$output, res231$output)
write.table(res.omegac2, "example results for omega c2-u1v05.txt")



res311 = repfun(N=300, n=23400, u=2, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res321 = repfun(N=300, n=23400, u=2,  v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res331 = repfun(N=300, n=23400, u=2, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac3 = rbind(res311$output, res321$output, res331$output)
write.table(res.omegac3, "example results for omega c3-u2v05.txt")


res411 = repfun(N=1000, n=23400, u=3, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res421 = repfun(N=1000, n=23400, u=3, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3, gam=0)
res431 = repfun(N=1000, n=23400, u=3, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3, gam=0)

res.omegac4 = rbind(res411$output, res421$output, res431$output)
write.table(res.omegac4, "example results for omega  c4-u3v05.txt")

res511 = repfun(N=300, n=23400, u=4, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res521 = repfun(N=300, n=23400, u=4,  v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res531 = repfun(N=300, n=23400, u=4, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac5 = rbind(res511$output, res521$output, res531$output)
write.table(res.omegac5, "example results for omega c5-u4v05.txt")

res611 = repfun(N=300, n=23400, u=0.5, v=1, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res621 = repfun(N=300, n=23400, u=0.5,  v=1, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res631 = repfun(N=300, n=23400, u=0.5, v=1, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac6 = rbind(res611$output, res621$output, res631$output)
write.table(res.omegac6, "example results for omega c6-u05v1.txt")

res711 = repfun(N=300, n=23400, u=1, v=1, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res721 = repfun(N=300, n=23400, u=1, v=1, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3, gam=0)
res731 = repfun(N=300, n=23400, u=1, v=1, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3, gam=0)

res.omegac7 = rbind(res711$output, res721$output, res731$output)
write.table(res.omegac7, "example results for omega  c7-u1v1.txt")

res811 = repfun(N=300, n=23400, u=2, v=1, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res821 = repfun(N=300, n=23400, u=2,  v=1, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res831 = repfun(N=300, n=23400, u=2, v=1, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac8 = rbind(res811$output, res821$output, res831$output)
write.table(res.omegac8, "example results for omega c8-u2v1.txt")

res911 = repfun(N=300, n=23400, u=3, v=1, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res921 = repfun(N=300, n=23400, u=3,  v=1, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res931 = repfun(N=300, n=23400, u=3, v=1, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac9 = rbind(res911$output, res921$output, res931$output)
write.table(res.omegac9, "example results for omega c9-u3v1.txt")


res1011 = repfun(N=300, n=23400, u=4, v=1, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res1021 = repfun(N=300, n=23400, u=4,  v=1, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res1031 = repfun(N=300, n=23400, u=4, v=1, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac10 = rbind(res1011$output, res1021$output, res1031$output)
write.table(res.omegac10, "example results for omega c10-u4v1.txt")


res1111 = repfun(N=300, n=23400, u=0.5, v=2, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res1121 = repfun(N=300, n=23400, u=0.5,  v=2, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res1131 = repfun(N=300, n=23400, u=0.5, v=2, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac11 = rbind(res1111$output, res1121$output, res1131$output)
write.table(res.omegac11, "example results for omega c11-u05v2.txt")


res1211 = repfun(N=300, n=23400, u=1, v=2, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res1221 = repfun(N=300, n=23400, u=1,  v=2, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res1231 = repfun(N=300, n=23400, u=1, v=2, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac12 = rbind(res1211$output, res1221$output, res1231$output)
write.table(res.omegac11, "example results for omega c12-u1v2.txt")

res1311 = repfun(N=300, n=23400, u=2, v=2, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res1321 = repfun(N=300, n=23400, u=2,  v=2, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res1331 = repfun(N=300, n=23400, u=2, v=2, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac13 = rbind(res1311$output, res1321$output, res1331$output)
write.table(res.omegac13, "example results for omega c13-u2v2.txt")


res1411 = repfun(N=300, n=23400, u=3, v=2, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res1421 = repfun(N=300, n=23400, u=3,  v=2, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res1431 = repfun(N=300, n=23400, u=3, v=2, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac14 = rbind(res1411$output, res1421$output, res1431$output)
write.table(res.omegac14, "example results for omega c14-u3v2.txt")


res1511 = repfun(N=300, n=23400, u=4, v=2, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res1521 = repfun(N=300, n=23400, u=4,  v=2, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res1531 = repfun(N=300, n=23400, u=4, v=2, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac15 = rbind(res1511$output, res1521$output, res1531$output)
write.table(res.omegac15, "example results for omega c15-u4v2.txt")


res1611 = repfun(N=300, n=23400, u=0.5, v=3, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res1621 = repfun(N=300, n=23400, u=0.5,  v=3, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res1631 = repfun(N=300, n=23400, u=0.5, v=3, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac16 = rbind(res1611$output, res1621$output, res1631$output)
write.table(res.omegac16, "example results for omega c16-u05v3.txt")



res1711 = repfun(N=300, n=23400, u=1, v=3, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res1721 = repfun(N=300, n=23400, u=1,  v=3, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res1731 = repfun(N=300, n=23400, u=1, v=3, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac17 = rbind(res1711$output, res1721$output, res1731$output)
write.table(res.omegac17, "example results for omega c17-u1v3.txt")



res1811 = repfun(N=300, n=23400, u=2, v=3, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res1821 = repfun(N=300, n=23400, u=2,  v=3, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res1831 = repfun(N=300, n=23400, u=2, v=3, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac18 = rbind(res1811$output, res1821$output, res1831$output)
write.table(res.omegac18, "example results for omega c18-u2v3.txt")

res1911 = repfun(N=300, n=23400, u=3, v=3, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res1921 = repfun(N=300, n=23400, u=3,  v=3, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res1931 = repfun(N=300, n=23400, u=3, v=3, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac19 = rbind(res1911$output, res1921$output, res1931$output)
write.table(res.omegac19, "example results for omega c19-u3v3.txt")



res2011 = repfun(N=300, n=23400, u=4, v=3, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res2021 = repfun(N=300, n=23400, u=4,  v=3, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res2031 = repfun(N=300, n=23400, u=4, v=3, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac20 = rbind(res2011$output, res2021$output, res2031$output)
write.table(res.omegac20, "example results for omega c20-u4v3.txt")


res2111 = repfun(N=300, n=23400, u=4, v=4, alpha1=0.5, alpha2=0.5, omega1=0.01, omega2=0.01, theta=1/3,   gam=0)
res2121 = repfun(N=300, n=23400, u=4,  v=4, alpha1=0.5, alpha2=0.5, omega1=0.03, omega2=0.03, theta=1/3,   gam=0)
res2131 = repfun(N=300, n=23400, u=4, v=4, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegac21 = rbind(res2111$output, res2121$output, res2131$output)
write.table(res.omegac21, "example results for omega c21-u4v4.txt")









########try omega=0.01  0.03  0.05


res131 = repfun(N=3000, n=23400, u=1, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.05, omega2=0.05, theta=1/3, gam=0)


res231 = repfun(N=3000, n=23400, u=1, v=1, alpha1=0.5, alpha2=0.5, omega1=0.05, omega2=0.05, theta=1/3,   gam=0)

res331 = repfun(N=3000, n=23400, u=1, v=1.5, alpha1=0.5, alpha2=0.5, omega1=0.05, omega2=0.05, theta=1/3,   gam=0)

res.omegar3 = cbind(res131$output, res231$output, res331$output)
write.table(res.omegar3, "example results for omega r3-omega005-4.txt")

res131 = repfun(N=3000, n=23400, u=1, v=0.5, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3, gam=0)


res231 = repfun(N=3000, n=23400, u=1, v=1, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res331 = repfun(N=3000, n=23400, u=1, v=1.5, alpha1=0.5, alpha2=0.5, omega1=0.1, omega2=0.1, theta=1/3,   gam=0)

res.omegar3 = cbind(res131$output, res231$output, res331$output)
write.table(res.omegar3, "example results for omega r3-omega01-3.txt")
