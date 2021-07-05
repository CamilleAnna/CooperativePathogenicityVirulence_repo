
mytheme<- theme_bw()+
  theme(legend.position="none",
        panel.border= element_blank(),
        axis.text.y = element_text(face="bold", colour="black", size=10),
        axis.text.x = element_text(face="bold", colour="black", size=11),
        axis.title.y = element_text(face="bold", colour="black", size=11),
        axis.title.x = element_text(face="bold", colour="black", size=11),
        axis.line.y = element_line(color="black", size = 0.5),
        axis.line.x = element_line(color="black", size = 0.5),
        plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))

mytheme_leg<- theme_bw()+
  theme(#legend.position="none",
    panel.border= element_blank(),
    axis.text.y = element_text(face="bold", colour="black", size=10),
    axis.text.x = element_text(face="bold", colour="black", size=11),
    axis.title.y = element_text(face="bold", colour="black", size=11),
    axis.title.x = element_text(face="bold", colour="black", size=11),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.line.x = element_line(color="black", size = 0.5),
    plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5))
