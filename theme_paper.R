
theme_paper <- function(){ 
  font <- "Arial"   #assign font family up front
  
  theme_classic2() %+replace%    #replace elements we want to change
    
    theme(
      text = element_text(colour = "black"
                          ),
      #grid elements
      panel.grid.major = element_blank(),  #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines       #strip axis ticks
      plot.background = element_blank(),
      #panel.border=element_rect(color="black", fill = NA),
      panel.background = element_blank(),
      
      #axis lines
      axis.line = element_line(colour = 'black', size = 0.5),
      axis.ticks = element_line(colour = "black", size = 0.5),
      
      
      
      #text elements
      
      axis.title = element_text(             #axis titles
        #family = font,            #font family
        size = 8,
        color = "black"),               #font size
      
      axis.text = element_text(              #axis text
        #family = font,            #axis famuly
        size = 7,
        color = "black"),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10),
        color = "black"),
      
      legend.text = element_text(              #axis text
        #family = font,            #axis famuly
        size = 7,
        color = "black"),  

      
      legend.title = element_blank(),
      legend.background = element_blank()

      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}
