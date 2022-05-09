# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 01:57:58 2021

@author: Kedar
"""

import pygame, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")

import matplotlib.backends.backend_agg as agg

import pointDefects as pd
import pylab

pygame.init()

X = 900  # screen width
Y = 600  # screen height

WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
RED = (255, 50, 50)
YELLOW = (255, 255, 0)
GREEN = (0, 255, 50)
BLUE = (50, 50, 255)
GREY = (200, 200, 200)
ORANGE = (200, 100, 50)
CYAN = (0, 255, 255)
MAGENTA = (255, 0, 255)
TRANS = (1, 1, 1)






# define class for sliders 

class Slider():
    def __init__(self, name, val, maxi, mini, pos):
        self.name = name
        self.val = val  # start value
        self.maxi = maxi  # maximum at slider position right
        self.mini = mini  # minimum at slider position left
        self.xpos = pos  # x-location on screen
        self.ypos = 525
        self.surf = pygame.surface.Surface((100, 50))
        self.hit = False  # the hit attribute indicates slider movement due to mouse interaction
        name = name + str(self.val)[0:5]
        self.txt_surf = font.render(name, 1, BLACK)
        self.txt_rect = self.txt_surf.get_rect(center=(50, 14))

        # Static graphics - slider background #
        self.surf.fill((100, 100, 100))
        pygame.draw.rect(self.surf, GREY, [0, 0, 100, 50], 3)
        pygame.draw.rect(self.surf, ORANGE, [10, 10, 80, 12], 0)
        pygame.draw.rect(self.surf, WHITE, [10, 30, 80, 5], 0)

        self.surf.blit(self.txt_surf, self.txt_rect)  # this surface never changes

        # dynamic graphics - button surface #
        self.button_surf = pygame.surface.Surface((20, 20))
        self.button_surf.fill(TRANS)
        self.button_surf.set_colorkey(TRANS)
        pygame.draw.circle(self.button_surf, BLACK, (10, 10), 6, 0)
        pygame.draw.circle(self.button_surf, ORANGE, (10, 10), 4, 0)

    
    def draw(self):
        """ Combination of static and dynamic graphics in a copy of
    the basic slide surface
    """
        # static
        surf = self.surf.copy()
        
        # dynamic

        pos = (10+int((self.val-self.mini)/(self.maxi-self.mini)*80), 33)
        self.button_rect = self.button_surf.get_rect(center=pos)
        surf.blit(self.button_surf, self.button_rect)
        self.button_rect.move_ip(self.xpos, self.ypos)  # move of button box to correct screen position
         
        
        # screen
        screen.blit(surf, (self.xpos, self.ypos),)
        #update text value of param
        pygame.draw.rect(self.surf, ORANGE, [10, 10, 80, 12], 0)
        self.txt_surf = font.render(self.name.split(":")[0]+": "+str(self.val)[0:5], 1, BLACK)    
        screen.blit(self.txt_surf, (self.xpos+20,self.ypos+8))#,self.txt_surf.get_rect(center=(50, 15)))

    def move(self):
        """
    The dynamic part; reacts to movement of the slider button.
    """
        # static
        
        self.val = (pygame.mouse.get_pos()[0] - self.xpos - 10) / 80 * (self.maxi - self.mini) + self.mini
        if self.val < self.mini:
            self.val = self.mini
        if self.val > self.maxi:
            self.val = self.maxi
        
        
font = pygame.font.SysFont("Verdana", 12)
screen = pygame.display.set_mode((X, Y))
clock = pygame.time.Clock()        




def main():
    
    # initialize sliders
    x = np.arange(8,900,(900)/8)
    G1 = Slider("G1: ", 56, 100, 0, x[0])
    G2 = Slider("G2: ", 260, 500, 0, x[1])
    tg = Slider("tg: ", 0, 2*np.pi, 0, x[2])
    pg = Slider("pg: ", 0, 2*np.pi, 0, x[3])
    te = Slider("te: ", 0, 2*np.pi, 0, x[4])
    pe = Slider("pe: ", 0, 2*np.pi, 0, x[5])
    tB = Slider("tB: ", 0.955, 2*np.pi, 0, x[6])
    f = Slider("f: ", 0, 10, 0, x[7])
    slides = [G1,G2,tg,pg,te,pe,tB,f]
    
    
    
    while True:
        #reset switch
        button = pygame.Rect(850, 440, 45, 45)
        text = font.render('RESET' , True , BLACK) 
        
        #draw slides
        for s in slides:
            s.draw()
        #Move slides
        for s in slides:
            if s.hit:
                s.move()
                
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif event.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                if button.collidepoint(pos):
                                    
                    main()
                    
                for s in slides:
                    if s.button_rect.collidepoint(pos):
                        s.hit = True
            elif event.type == pygame.MOUSEBUTTONUP:
                for s in slides:
                    s.hit = False
            
        #draw plots to pygame through bitmap
        fig = pylab.figure(num=1,figsize=[9, 5], # Inches
                           dpi=100, clear=True        # 100 dots per inch, so the resulting buffer is 400x400 pixels
                           )
        B, EEEvg, EEEve, Eg, Ee = pd.SiVModel(G1.val,G2.val,tg.val,pg.val,te.val,pe.val,tB.val,f.val)
        ax = fig.gca()
        
        ax.plot(B,EEEvg[:][:] , label = "Ground States")
        ax.plot(B,EEEve[:][:] , label = "Excited States")
        
        pylab.xlabel('magnetic field (T)')
        pylab.ylabel('Transition Frequency (Hz)')
        pylab.title('Calculated splitting of SiV electronic levels for increasing magnetic field')
        pylab.legend()
        canvas = agg.FigureCanvasAgg(fig)
        
        canvas.draw()
        renderer = canvas.get_renderer()
        raw_data = renderer.tostring_rgb()
        screen = pygame.display.get_surface()
    
        size = canvas.get_width_height()
        #draw plot
        surf = pygame.image.fromstring(raw_data, size, "RGB")
        screen.blit(surf, (0,0))
        #draw reset button
        pygame.draw.rect(screen,[255, 0, 0], button,)
        screen.blit(text , (button.x+1/5*button.width/2,button.y+1/2*button.height/2))#  draw button
        #update screen
        pygame.display.flip()
        #setp animation timestep
        clock.tick(60) #fps>60 
main()