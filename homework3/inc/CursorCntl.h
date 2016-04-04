#ifndef CURSOR_CONTROL
#define CURSOR_CONTROL

//#define DEBUG

typedef enum {
	BLACK,
	BOLDBLACK,
	RED,
	ORANGE,
	GREEN,
	BOLDGRAY,
	YELLOW,
	BOLDGRAY2,
	BLUE,
	BOLDSTANDARD,
	MAGENTA,
	PURPLE,
	CYAN,
	WHITE,
	BOLDWHITE,
	RESET
}Colors;

void setColor( Colors color);

void gotoXY(int X, int Y);

void saveCursor();

void clearEOL();

void clearPAGE();

void recallCursor();

void clearBelowLine(int line);

void clearLine(int line);

void printError(char * message);

void printOK(char * message);

void printWarning(char * message);

void printPrompt(char * message);


#endif
