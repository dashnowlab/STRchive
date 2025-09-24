import { cloneElement, Fragment, useMemo } from "react";
import { FaArrowDown, FaArrowUp, FaPlus, FaTrash } from "react-icons/fa6";
import { compileSchema, draft2020, extendDraft } from "json-schema-library";
import {
  cloneDeep,
  isEmpty,
  isObject,
  mapValues,
  range,
  uniq,
} from "lodash-es";
import { get, join, remove, set, split } from "@sagold/json-pointer";
import Button from "@/components/Button";
import ComboBox from "@/components/ComboBox";
import Heading from "@/components/Heading";
import Help from "@/components/Help";
import NumberBox from "@/components/NumberBox";
import Select from "@/components/Select";
import TextBox from "@/components/TextBox";
import { makeList } from "@/util/format";
import classes from "./SchemaForm.module.css";

/** form automatically generated from json schema */
const SchemaForm = ({ schema, sections, data, onChange, children }) => {
  /** compile schema */
  const rootNode = useMemo(
    () => compileSchema(schema, { drafts: [draft] }),
    [schema],
  );

  /** if no initial data, generate blank values */
  data ??= rootNode.getData(null, { extendDefaults: false });

  /** revalidate when data changes */
  const { errors } = useMemo(() => {
    rootNode.context.data = data;
    return rootNode.validate(data);
  }, [rootNode, data]);

  return (
    <>
      {children && <section>{children}</section>}

      {sections.map((section, index) => (
        <section key={index}>
          <Heading level={2}>{section}</Heading>

          <Field
            rootNode={rootNode}
            schema={schema}
            section={section}
            node={rootNode}
            data={data}
            setData={onChange}
            errors={errors}
          />
        </section>
      ))}
    </>
  );
};

export default SchemaForm;

/** add custom keyword for validating one value less than other value */
const lessThan = {
  id: "lessThan",
  keyword: "less_than",
  addValidate: () => true,
  validate: ({ node, data, pointer }) => {
    const thisValue = data;
    const fullData = node.context.data;
    const otherKey = node.schema.less_than;
    const otherValue = get(fullData, otherKey);
    if (thisValue === null || otherValue === null) return;
    if (thisValue <= otherValue) return;
    return node.createError("compare-value-error", {
      pointer,
      value: otherValue,
      schema: { ...node.schema, message: `Must be less than ${otherKey}` },
    });
  },
};

/** add custom keyword for validating one value greater than other value */
const greaterThan = {
  id: "greaterThan",
  keyword: "greater_than",
  addValidate: () => true,
  validate: ({ node, data, pointer }) => {
    const thisValue = data;
    const fullData = node.context.data;
    const otherKey = node.schema.greater_than;
    const otherValue = get(fullData, otherKey);
    if (thisValue === null || otherValue === null) return;
    if (thisValue >= otherValue) return;
    return node.createError("compare-value-error", {
      pointer,
      value: otherValue,
      schema: { ...node.schema, message: `Must be greater than ${otherKey}` },
    });
  },
};

/** https://www.npmjs.com/package/json-schema-library#user-content-extending-a-draft */
const draft = extendDraft(draft2020, {
  /** match all schema */
  $schemaRegEx: ".*",
  /** add custom validation rules */
  keywords: [lessThan, greaterThan],
});

/**
 * @param {Object} props
 * @param {import("json-schema-library").SchemaNode} props.node
 * @param {import("json-schema-library").JsonError[]} props.errors
 */
const Field = ({
  rootNode,
  schema,
  section,
  node,
  path = "",
  data,
  setData,
  errors,
}) => {
  /** are we at top level of schema */
  const level = split(path).length;

  /** filter out fields that should not be displayed in this section */
  if (level === 1 && node.schema.section !== section) return;

  /** schema props */
  const {
    title,
    description,
    placeholder,
    examples,
    type,
    pattern,
    enum: _enum,
    enum_descriptions,
    minimum,
    maximum,
    less_than,
    greater_than,
    hide,
    multiline,
    combobox,
  } = node.schema;

  /** explicitly hide field */
  if (hide) return;

  /** normalize type to array */
  const types = [type].flat();

  /** form data name */
  const name = path;

  /** get nested data value from path */
  const value = get(data, path);

  /** set nested data value from path */
  const onChange = (value) => {
    let _data = cloneDeep(data);
    _data = set(_data, path, value);
    setData(_data);
  };

  /** help tooltip to show */
  const tooltip = [
    description,
    ...(examples?.length
      ? [
          "Examples:",
          `<ul>${examples.map((ex) => `<li>${ex}</li>`).join("")}</ul>`,
        ]
      : []),
    ...(!isEmpty(enum_descriptions) && value in enum_descriptions
      ? [`${value}: ${enum_descriptions[value]}`]
      : []),
  ]
    .filter(Boolean)
    .map((line) => `<div>${line}</div>`)
    .join("");

  /** is the field required to be set */
  const required =
    schema.required?.includes(split(name).pop()) && !types.includes("null");

  /** full label elements to show */
  const label = (
    <>
      {title ?? name ?? ""}
      {tooltip && <Help>{tooltip}</Help>}
    </>
  );

  /** get validation errors associated with this field */
  const fieldErrors = errors.filter(
    ({ data }) => data.pointer === path.replace(/^#?\/?/, "#/"),
  );
  /** is there an error on this field */
  const isError = !!fieldErrors.length;

  /** error component to show */
  const error = isError ? (
    <span className={classes.error}>
      {fieldErrors
        .map(({ code, data, message }) => {
          console.log(code, data, message);
          /** prettify error message */
          if (code === "pattern-error")
            return (
              <>
                Must match pattern <code>{data.pattern}</code>
              </>
            );
          if (code === "enum-error")
            return <>Must be {makeList(data.schema.enum, "code")}</>;
          if (code === "compare-value-error") return <>{data.schema.message}</>;
          if (code === "type-error") {
            if (
              [data.expected].flat().includes("string") &&
              (data.value === "" || data.value == null)
            )
              return <>Must not be blank</>;
            else return <>Must be {makeList(data.expected, "code")}</>;
          }
          if (code === "maximum-error")
            return (
              <>
                Maximum <code>{data.maximum}</code>
              </>
            );
          if (code === "minimum-error")
            return (
              <>
                Minimum <code>{data.minimum}</code>
              </>
            );
          if (code === "unique-items-error") return <>Must not be duplicate</>;
          if (code === "min-items-error")
            return (
              <>
                Must have at least <code>{data.minItems}</code> items
              </>
            );
          if (code === "max-items-error")
            return (
              <>
                Must have at most <code>{data.maxItems}</code> items
              </>
            );
          return message;
        })
        .map((error, index) => (
          <Fragment key={index}>{error}. </Fragment>
        ))}
    </span>
  ) : null;

  /** field control to show */
  let control = <span className={classes.error}>missing control</span>;

  /** ref function */
  const ref = (el) => {
    if (!el || !("setCustomValidity" in el)) return;
    /** also set native browser form validity to get various built-in features */
    if (isError) el?.setCustomValidity("Error with this field");
    else el?.setCustomValidity("");
  };

  if (types.includes("object"))
    /** object group */
    control = (
      <div className={classes.object}>
        {level > 0 && <div className={classes.heading}>{label}</div>}
        {Object.keys(node.schema.properties).map((key) => {
          return (
            <Field
              key={key}
              rootNode={rootNode}
              schema={schema}
              section={section}
              node={node.getChildSelection(key)[0]}
              path={`${path}/${key}`}
              data={data}
              setData={setData}
              errors={errors}
            />
          );
        })}
      </div>
    );
  else if (types.includes("array")) {
    /** array group */
    const items = get(data, path)?.length ?? 0;

    control = (
      <div className={classes.array}>
        <div className={classes.heading}>{label}</div>
        {range(items).map((index) => (
          <Fragment key={index}>
            <Field
              key={index}
              schema={schema}
              section={section}
              node={node.getChildSelection()[0]}
              path={`${path}/${index}`}
              data={data}
              setData={setData}
              errors={errors}
            />
            <div className={classes.actions}>
              <Button
                disabled={index === 0}
                data-tooltip="Move up"
                onClick={() => {
                  let _data = cloneDeep(data);
                  const aPath = join(path, String(index - 1));
                  const bPath = join(path, String(index));
                  const aValue = get(data, aPath);
                  const bValue = get(data, bPath);
                  _data = set(_data, aPath, bValue);
                  _data = set(_data, bPath, aValue);
                  setData(_data);
                }}
              >
                <FaArrowUp />
              </Button>
              <Button
                disabled={index === items - 1}
                data-tooltip="Move down"
                onClick={() => {
                  let _data = cloneDeep(data);
                  const aPath = join(path, String(index));
                  const bPath = join(path, String(index + 1));
                  const aValue = get(data, aPath);
                  const bValue = get(data, bPath);
                  _data = set(_data, aPath, bValue);
                  _data = set(_data, bPath, aValue);
                  setData(_data);
                }}
              >
                <FaArrowDown />
              </Button>
              <Button
                data-tooltip="Remove"
                onClick={() => {
                  let _data = cloneDeep(data);
                  _data = remove(_data, join(path, String(index)));
                  setData(_data);
                }}
              >
                <FaTrash />
              </Button>
            </div>
          </Fragment>
        ))}

        <Button
          style={{ gridColumn: "1 / -1" }}
          design="plain"
          onClick={() => {
            /** add new item to array */
            let _data = cloneDeep(data);
            /** get child schema */
            let newItem = node
              .getNodeChild("items")
              .node.getData(undefined, {});
            /** modify default value fills */
            if (newItem === "") newItem = null;
            if (isObject(newItem)) newItem = mapValues(newItem, () => null);
            /** insert item */
            _data = set(_data, [...split(path), "[]"], newItem);
            setData(_data);
          }}
        >
          <FaPlus />
          <span>Add</span>
        </Button>
      </div>
    );
  } else if (_enum?.length) {
    /** single select */
    const options = uniq(["", ..._enum.filter(Boolean)]).map((value) => ({
      value,
      label: value,
      tooltip: enum_descriptions?.[value],
    }));
    control = (
      <Select
        options={options}
        value={value || ""}
        onChange={(value) => onChange(value || null)}
      />
    );
  } else if (types.includes("string"))
    if (combobox && examples) {
      /** text input with autocomplete suggestions */
      control = (
        <ComboBox
          options={examples.map((example) => ({
            value: example,
            label: example,
          }))}
          value={value}
          onChange={onChange}
        />
      );
    } else
      /** plain text input */
      control = (
        <TextBox
          placeholder={
            placeholder ?? (examples?.length ? `Ex: ${examples[0]}` : "")
          }
          multi={multiline}
          // pattern={pattern}
          minLength={required ? 0 : undefined}
          value={value || ""}
          onChange={(value) => onChange(value || null)}
        />
      );
  else if (types.includes("number") || types.includes("integer")) {
    /** min value allowed to be input */
    const min = Math.max(
      ...[minimum, get(data, greater_than)].filter(
        (value) => typeof value === "number",
      ),
    );
    /** max value allowed to be input */
    const max = Math.min(
      ...[maximum, get(data, less_than)].filter(
        (value) => typeof value === "number",
      ),
    );
    control = (
      <NumberBox
        step={types.includes("integer") ? 1 : "any"}
        min={Number.isFinite(min) ? min : undefined}
        max={Number.isFinite(max) ? max : undefined}
        value={value}
        onChange={onChange}
      />
    );
  }

  return (
    <Fragment>
      {typeof control.type === "function"
        ? cloneElement(
            control,
            /** props that should be on all component controls */
            { ref, name, label, required },
          )
        : control}
      {error}
    </Fragment>
  );
};
