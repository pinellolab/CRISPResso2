'''
This file is derived from https://github.com/mikeckennedy/jinja_partials and is subject to the following license:

MIT License

Copyright (c) 2021 Michael Kennedy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

from functools import partial

from markupsafe import Markup

try:
    import flask
except ImportError:
    flask = None


def render_partial(template_name, renderer=None, markup=True, **data):
    """
    Renders a partial template and returns the result. If `markup` is True, the result is wrapped in a `Markup` object.
    """
    if renderer is None:
        if flask is None:
            raise PartialsException('No renderer specified')
        renderer = flask.render_template

    if markup:
        return Markup(renderer(template_name, **data))

    return renderer(template_name, **data)


def generate_render_partial(renderer, markup=True):
    """
    Returns a partial function that renders a template using the specified renderer.
    """
    return partial(render_partial, renderer=renderer, markup=markup)
