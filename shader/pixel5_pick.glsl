
uniform vec4 ObjectId;
out vec4 color;

// layout(early_fragment_tests) in;
// ^- apparently not supported in OpenGL ES.

void main(void)
{
   color = ObjectId;
}
