in vec4 in_Position;

out vec4 PixelPos;

void main() {
    PixelPos = in_Position;
    gl_Position = in_Position;
}
