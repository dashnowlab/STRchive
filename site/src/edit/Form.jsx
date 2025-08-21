import { Fragment, useCallback, useState } from "react";
import { LuSend } from "react-icons/lu";
import { groupBy, upperFirst } from "lodash-es";
import Button from "@/components/Button";
import FormWrapper from "@/components/Form";
import Heading from "@/components/Heading";
import Help from "@/components/Help";
import NumberBox from "@/components/NumberBox";
import Select from "@/components/Select";
import TextBox from "@/components/TextBox";
import { validate } from "@/data/schema";
import classes from "./Form.module.css";
import {
  properties,
  required as requiredKeys,
} from "~/STRchive-loci.schema.json";

/** edit/new locus page */

/** split fields in schema by manually labeled "section" */
const sections = groupBy(
  Object.entries(properties).map(([key, value]) => ({ key, ...value })),
  "section",
);

const Form = ({ data: initialData = {} }) => {
  /** current form data */
  const [data, _setData] = useState(initialData);
  /** current error data */
  const [errors, _setErrors] = useState(null);

  /** set error func */
  const setErrors = useCallback(
    (newData) => {
      let errors = validate(newData ?? data) ?? [];
      errors = errors.map(({ instancePath, message, params }) => [
        instancePath.split("/").pop() || params.missingProperty,
        message,
      ]);
      errors = Object.fromEntries(errors);
      _setErrors(errors);
    },
    [data],
  );

  /** set data func */
  const setData = useCallback((key, value) => {
    _setData((data) => {
      const newData = { ...data, [key]: value };
      setErrors(newData);
      return newData;
    });
  }, []);

  return (
    <FormWrapper onSubmit={console.info}>
      <div></div>
      {Object.entries(sections).map(([section, properties]) => {
        if (section === "undefined") return <Fragment key={section} />;
        return (
          <section key={section}>
            <Heading level={2}>{section}</Heading>

            <div className={classes.form}>
              {properties.map((property) => {
                /** get different fields on property */
                const {
                  key,
                  type,
                  title,
                  description,
                  enum: _enum,
                  examples,
                  pattern,
                  minimum,
                  maximum,
                  less_than,
                  greater_than,
                } = property;

                /** get error message associated with this property, if any */
                const errorMessage = errors?.[key] || "";
                /** is there an error on this property */
                const error = !!errorMessage;
                /** error component to show */
                const errorRow = error ? (
                  <span className={classes.error}>{errorMessage}</span>
                ) : null;

                /** ref function */
                const ref = (el) => {
                  /** even though we show our own error messages, also set
                   * native browser form validity */
                  if (error) el?.setCustomValidity(errorMessage);
                  else el?.setCustomValidity("");
                };

                /** normalize type to array */
                const types = [type].flat();

                /** is null an allowable value */
                const allowNull = types.includes("null");

                /** is the field required to be set */
                const required = requiredKeys.includes(key) && !allowNull;

                /** help tooltip to show */
                const tooltip = [
                  description,
                  ...(examples?.length ? [" "] : []),
                  examples?.length ? `Examples: ${examples.join(", ")}` : "",
                ]
                  .filter(Boolean)
                  .join("<br/>");

                /** label component to show */
                const label = (
                  <>
                    {title ?? key}
                    {description && <Help>{tooltip}</Help>}
                  </>
                );

                /** current value of property in form data */
                const value = data[key] ?? "";

                /** on change handler for plain values */
                const onChange = (value) => {
                  if (value === "" && allowNull) setData(key, null);
                  else if (value === "" && !allowNull) setData(key, undefined);
                  else setData(key, value);
                };

                /** final field component to show */
                let field = <></>;

                if (_enum)
                  /** single select */
                  field = (
                    <Select
                      key={key}
                      ref={ref}
                      name={key}
                      label={label}
                      required={required}
                      options={_enum.map((value) => ({
                        value,
                        label: upperFirst(value),
                      }))}
                      value={value}
                      onChange={onChange}
                    />
                  );
                else if (types.includes("string"))
                  /** text box */
                  field = (
                    <TextBox
                      key={key}
                      ref={ref}
                      name={key}
                      label={label}
                      required={required}
                      pattern={pattern}
                      value={value}
                      onChange={onChange}
                    />
                  );
                else if (
                  types.includes("integer") ||
                  types.includes("number")
                ) {
                  /** number box */

                  /** min value allowed to be input */
                  const min = Math.max(
                    ...[minimum, data[greater_than]].filter(
                      (value) => typeof value === "number",
                    ),
                  );
                  /** max value allowed to be input */
                  const max = Math.min(
                    ...[maximum, data[less_than]].filter(
                      (value) => typeof value === "number",
                    ),
                  );
                  field = (
                    <NumberBox
                      key={key}
                      ref={ref}
                      name={key}
                      label={label}
                      required={required}
                      pattern={pattern}
                      step={types.includes("integer") ? 1 : undefined}
                      min={min}
                      max={max}
                      value={value}
                      onChange={onChange}
                    />
                  );
                } else return <Fragment key={key} />;

                return (
                  <Fragment key={key}>
                    {field}
                    {errorRow}
                  </Fragment>
                );
              })}
            </div>
          </section>
        );
      })}

      <section>
        <Button type="submit" design="bubble" onClick={() => setErrors()}>
          <LuSend />
          <span>Submit</span>
        </Button>
      </section>
    </FormWrapper>
  );
};

export default Form;
