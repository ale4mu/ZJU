#include "adventure.h"

int main(){
    Game game;
    cout << "Do u want to fix the map?(y/n)" << endl;
    if(getchar() == 'y')
        game.set_fixed_map();
    game.set_game();
    game.play();
}